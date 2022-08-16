#include <vector>
#include <random>
#include <tuple>
#include "timer.cpp"
#include "omp.h"
namespace dp {

    using Vector = std::vector<double>;
    using Matrix = std::vector<Vector>;

    // initialization utilities

    Vector full(int n, double fill_value){
        Vector x(n);
        fill(x.begin(), x.end(), fill_value);
        return x;
    }

    Matrix full(int n, int m, double fill_value) {
        Matrix x(m);
        fill(x.begin(), x.end(), full(n, fill_value));
        return x;
    }

    auto ones(int n) {
        return full(n, 1);
    }

    auto ones(int n, int m) {
        return full(n,m,1);
    }

    auto zeros(int n) {
        return full(n,0);
    }

    auto zeros(int n, int m) {
        return full(n,m,0);
    }

    auto rand(int n) {
        Vector x;
        for (int i = 0; i < n; i++) {
            x.push_back(double(std::rand() / double(1)));
        }
        return x;
    }

    auto rand(int n, int m) {
        Matrix x(m);
        std::fill(x.begin(), x.end(), rand(n));
        return x;
    }

    Matrix diagonal(const int n, const double fill_value) {
        auto x = zeros(n,n);

        for (int i = 0; i < n; i++)
            x[i][i] = fill_value;
        return x;
    }

    // vector operations

    double dot(Vector x, Vector y, unsigned long long start=0, unsigned long long end=-1) {
        double product = 0;

        if (end==-1)
            end = x.size();

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

        for(int i=0; i<A.size();i++)
            for(int j=0; j<A[0].size(); j++)
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
        Matrix A = sum(ones(n, n),diagonal(n, n), nw);
        Vector b = full(n,2*n);
        return std::make_tuple(A, b);
    }
}
