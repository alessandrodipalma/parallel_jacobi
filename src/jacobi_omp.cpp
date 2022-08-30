#include "../include/solvers.h"
#include <omp.h>
#include <numeric>

Vector jacobi_omp::solve(Matrix A, Vector b, const int max_iter, int nw,
                         const std::function<bool(Vector &)> stopping_criteria) {

    int n = A.size();

    Vector diag(n);
    for (int i = 0; i < n; i++) {
        diag[i] = A[i][i];
        A[i][i] = 0;
    }

    Vector x(n, 0);
    Vector x_new(n, 0);

    int k = 0;

    while (k < max_iter) {
        #pragma omp parallel for num_threads(nw)
        for (int i = 0; i < n; i++) {
            double s = std::inner_product(A[i].begin(), A[i].end(), x.begin(), 0.0);
            x_new[i] = (b[i] - s) / diag[i];
        }

        if (stopping_criteria(x_new)) {
            std::cout << "in " << k << " iter" << std::endl;
            k=max_iter;
        } else {
            x.swap(x_new);
            k++;
        }
    }
    return x_new;
}
