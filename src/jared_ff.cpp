#include "../include/solvers.h"
#include "../include/timer.hpp"
#include <ff/parallel_for.hpp>

#include <thread>
#include <numeric>

Vector jared_ff::solve(Matrix A, const Vector b, const int max_iter, int nw,
                              std::function<bool(Vector &)> stopping_criteria) {

    int n = A.size();
    Vector diag(n);
    for (int i = 0; i < n; i++) {
        diag[i] = A[i][i];
        A[i][i] = 0;
    }

    Vector x(n, 0);
    Vector x_new(n, 0);
    Matrix sigma_parts = zeros(n, nw);
    Vector sigma(n, 0);
    const Vector zero(n,0);

    auto body = [&A, &diag, &b, &x, &x_new, &sigma_parts, n](const long i, int tid) {

        double s = std::inner_product(A[i].begin(), A[i].end(), x.begin(), 0.0);
        x_new[i] = (b[i] - s) / diag[i];

        for (int j = 0; j < n; j++) {
            sigma_parts[tid][j] += x_new[i] * A[i][j];
        }
    };

    auto f = [&](const long j, Vector &sum) {
        sum = sum + sigma_parts[j];
        std::fill(sigma_parts[j].begin(), sigma_parts[j].end(), 0.0);
    };
    auto fred = [](Vector &sum, const Vector partial_sum) {
        sum = sum + partial_sum;
    };


    ff::ParallelForReduce<Vector> pf(nw, false, true);
    for (int k = 0; k < max_iter; k += 2) {
        #ifdef MEASURE_ITERATES
        timer t("jared_ff iterate");
        #endif
        {
            #ifdef MEASURE_ITERATES
            timer t("jared_ff inner prducts");
            #endif
            pf.parallel_for_thid(0, n, 1, 1, body, nw);
        }


        if (stopping_criteria(x_new)) {
            #ifdef PRINT_ITER
            std::cout << "in " << k << " iter" << std::endl;
            #endif
            k = max_iter;
        } else {
            std::fill(sigma.begin(), sigma.end(), 0.0);
            {
                #ifdef MEASURE_ITERATES
                timer t("jared_ff reduce");
                #endif
                pf.parallel_reduce(sigma, zero, 0, nw, f, fred, nw);
            }


            x_new = (b - sigma) / diag;

            if (stopping_criteria != nullptr && stopping_criteria(x_new)) {
                #ifdef PRINT_ITER
                std::cout << "in " << k << " iter" << std::endl;
                #endif
                k = max_iter;
            } else {
                x.swap(x_new);
            }
        }


    }
    return x_new;
}

