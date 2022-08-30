#include <iostream>
#include "../include/solvers.h"

Vector jacobi_seq_separate_iter(Matrix A, Vector b, int max_iter, std::function<bool(Vector &)> stopping_criteria) {

    Vector x(b.size(),0);
    Vector x_new(b.size(),0);
    int n= A.size();
    Vector diag(n);
    for (int i = 0; i<n; i++) {
        diag[i] = A[i][i];
        A[i][i]=0;
    }

    for (int k = 0; k < max_iter; k++) {
        for (int i = 0; i < n; i++) {
            double s = std::inner_product(A[i].begin(), A[i].end(), x.begin(),0.0);
            x_new[i] = (b[i] - s) / diag[i];
        }
        if (stopping_criteria(x_new)) {
            std::cout << "in " << k << " iter" << std::endl;
            k=max_iter;
        }
        x.swap(x_new);
    }
    return x_new;
}
