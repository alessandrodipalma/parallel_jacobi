#include "../include/solvers.h"
#include <ff/parallel_for.hpp>
#include <thread>
#include <numeric>

Vector jacobi_ff::solve(Matrix A, const Vector b, const int max_iter, int nw,
                        std::function<bool(Vector &)> stopping_criteria) {

    int n = A.size();
    Vector diag(n);
    for (int i = 0; i < n; i++) {
        diag[i] = A[i][i];
        A[i][i] = 0;
    }

    Vector x(n, 0);
    Vector x_new(n, 0);

    int k = 0;

    ff::ParallelFor parallelFor(nw);

    while (k < max_iter) {
        parallelFor.parallel_for(0, n, 1, [&A, &diag, &b, &x, &x_new, k](const long i) {
                double s = std::inner_product(A[i].begin(), A[i].end(), x.begin(), 0.0);
                x_new[i] = (b[i] - s) / diag[i];
        }, nw);

        if (stopping_criteria != nullptr && stopping_criteria(x_new)) {
#ifdef PRINT_ITER
            std::cout << "in " << k << " iter" << std::endl;
#endif
            k = max_iter;
            return x_new;
        } else {
            x.swap(x_new);
            k++;
        }

    }
    return x_new;
}

