#include "../include/solvers.h"
#include <omp.h>
#include <numeric>

Vector jacobi_alternative::solve(Matrix A, Vector b, const int max_iter, int nw,
                         const std::function<bool(Vector &)> stopping_criteria) {

    int n = A.size();

    Vector diag(n);
    for (int i = 0; i < n; i++) {
        diag[i] = A[i][i];
        A[i][i] = 0;
    }


    Matrix sigma_parts = zeros(n, nw);
    Vector sigma(n,0);
    const Vector zero(nw, 0);
    Vector x(n, 0);
    Vector x_new(n, 0);

#pragma omp declare reduction(vector_plus : Vector : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

    for(int k = 0; k<max_iter; k+=2) {
        #pragma omp parallel for num_threads(nw)
        for (int i = 0; i < n; i++) {
            double s = std::inner_product(A[i].begin(), A[i].end(), x.begin(), 0.0);
            x_new[i] = (b[i] - s) / diag[i];

            int tid = omp_get_thread_num();
            for(int j=0; j<n; j++) {
                sigma_parts[tid][j] += x_new[i] * A[i][j];
            }
        }

        if (stopping_criteria(x_new)) {
#ifdef PRINT_ITER
            std::cout << "in " << k << " iter" << std::endl;
#endif
            k=max_iter;
        } else {
            std::fill(sigma.begin(), sigma.end(), 0.0);
#pragma omp parallel for shared(sigma_parts) reduction(vector_plus:sigma) num_threads(nw)
            for (int j = 0; j < nw; j++) {
                sigma = sigma + sigma_parts[j];
                std::fill(sigma_parts[j].begin(), sigma_parts[j].end(), 0.0);
            }

            x_new = (b - sigma) / diag;
        }

        if (stopping_criteria(x_new)) {
#ifdef PRINT_ITER
            std::cout << "in " << k << " iter" << std::endl;
#endif
            k=max_iter;
        } else {
            x.swap(x_new);
        }
    }
    return x_new;
}
