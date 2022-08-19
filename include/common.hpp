#ifndef PARALLEL_JACOBI_COMMON_H
#define PARALLEL_JACOBI_COMMON_H

#include <vector>
#include <random>
#include <tuple>
#include "omp.h"
#include "timer.hpp"

namespace dp {

    using Vector = std::vector<double>;
    using Matrix = std::vector<Vector>;

    // initialization utilities
    Matrix full(int n, int m, double fill_value) {
        Matrix x(m,Vector(n, fill_value));
        return x;
    }

    Matrix ones(int n, int m) {
        return full(n,m,1);
    }

    Matrix zeros(int n, int m) {
        return full(n,m,0);
    }

    Vector rand(int n) {
        Vector x;
        for (int i = 0; i < n; i++) {
            x.push_back(double(std::rand() / double(1)));
        }
        return x;
    }

    Matrix rand(int n, int m) {
        Matrix x(m, rand(n));
        return x;
    }

    Matrix diagonal(const int n, const double fill_value, int nw = 1) {
        auto x = zeros(n,n);
        #pragma omp parallel for num_threads(nw)
        for (int i = 0; i < n; i++)
            x[i][i] = fill_value;
        return x;
    }

    // vector operations
    double dot(Vector x, Vector y, unsigned long long start=0, unsigned long long end=-1, int nw=1) {
        double product = 0;

        if (end==-1)
            end = x.size();
        #pragma omp parallel for num_threads(nw)
        for (unsigned long long  i = start; i < end; i++) {
            product += x[i] * y[i];
        }
        return product;
    }

    Vector dot(Matrix A, Vector x) {
        Vector product(A.size());

        for (int i = 0; i < A.size(); i++) {
            product.assign(i, dot(A.at(i), x));
        }
        return product;
    }

    Matrix sum(const Matrix A, const Matrix B, const int nw=1) {
        auto n = A.size();
        auto m = A[0].size();
        Matrix res = zeros(n,m);

        for(int i=0; i<n;i++)
            for(int j=0; j<m; j++)
                res[i][j] = A[i][j] + B[i][j];
        return res;
    }

    // OPERATORS OVERLOADING
    Matrix operator+(const Matrix& A,const Matrix& B) {return sum(A,B);}
    double operator*(const Vector& x,const Vector& y) {return dot(x,y);}
    Vector operator*(const Matrix& A,const Vector& x) {return dot(A,x);}

    // checks
    bool are_ones(const Vector x) {
        for (auto &it: x) {
            if(it!=1){
                return false;
            }
        };
        return true;
    }

    auto generate_diagonally_dominant_problem(int n, int nw=1) {
        Matrix A = sum(ones(n, n),diagonal(n, n, nw), nw);
        Vector b(n,2*n);
        return std::make_tuple(A, b);
    }
}

#endif //PARALLEL_JACOBI_COMMON_H