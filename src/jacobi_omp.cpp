#include "../include/solvers.h"
#include <omp.h>
#include <numeric>

dp::Vector jacobi_omp(Matrix A, const Vector b, const int max_iter, int nw,
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

    #pragma omp parallel num_threads(nw)
    while (k < max_iter) {
        #pragma omp for
        for (int i = 0; i < n; i++) {
            double s = std::inner_product(A[i].begin(), A[i].end(), x.begin(), 0.0);
            x_new[i] = (b[i] - s) / diag[i];
        }

        #pragma omp critical
        {
            if (stopping_criteria(x_new)) {
                std::cout << "in " << k << " iter" << std::endl;
                return x_new;
            }
            x = x_new;
            k++;
        };

    }
    return x;
}
