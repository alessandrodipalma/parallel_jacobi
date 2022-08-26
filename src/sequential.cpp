#include <iostream>
#include "../include/solvers.h"

Vector jacobi_seq(Matrix A, Vector& b, int max_iter, int nw= 1) {

    Vector x(b.size(),0);

    for (int k = 0; k < max_iter; k++) {
        for (int i = 0; i < A.size(); i++) {
            double s = 0;

            for (int j = 0; j < A.size(); j++) {
                if (j != i)
                    s += A[i][j] * x[j];
            }
            x[i] = (b[i] - s) / A[i][i];
        }
        if (are_ones(x)){
            std::cout << "computed in " << k << " iterations" << std::endl;
            return x;
        }
    }
    return x;
}



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
            return x_new;
        }
        x = x_new;
    }
    return x;
}


Vector jacobi_prod(Matrix A, Vector b, int max_iter, int nw=1) {
    Vector x(b.size(), 0);
//    Vector x_new(b.size());
    for (int k = 0; k < max_iter; k++) {
        for (int i = 0; i < A.size(); i++) {
            double s1,s2;
            #pragma omp sections
            {
                #pragma omp section
                {
                    s1 = dot(A[i], x, 0, i, nw / 2);
                }
                #pragma omp section
                {
                    s2 = dot(A[i], x, i + 1, x.size(), nw / 2);
                }
            }

            x[i] = (b[i] - s1-s2) / A[i][i];
        }
    }
    return x;

}
