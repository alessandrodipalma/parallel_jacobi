#include "../../include/solvers.hpp"

using namespace dp;

Vector jacobi_seq(Matrix A, Vector b, int max_iter, int nw= 1) {

    Vector x = dp::rand(b.size());

    for (int k = 0; k < max_iter; k++) {
        for (int i = 0; i < A.size(); i++) {
            double s = 0;

            for (int j = 0; j < A.size(); j++) {
                if (j != i)
                    s = s + A[i][j] * x[j];
            }
            x[i] = (b[i] - s) / A[i][i];
        }
//        x = x_new;
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
                #pragma omp section/
                {
                    s2 = dot(A[i], x, i + 1, x.size(), nw / 2);
                }
            }

            x[i] = (b[i] - s1-s2) / A[i][i];
        }
    }
    return x;

}


