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

    for(int k = 0; k<max_iter; k++) {
        #pragma omp parallel for num_threads(nw)
        for (int i = 0; i < n; i++) {
            double s = std::inner_product(A[i].begin(), A[i].end(), x.begin(), 0.0);
            x_new[i] = (b[i] - s) / diag[i];
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
