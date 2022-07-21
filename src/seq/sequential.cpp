#include <iostream>
#include "../common/common.hpp"

using namespace dp;


Vector jacobi(Matrix A, Vector b, int max_iter) {

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
    }
    return x;
}


int main() {
    auto tup = generate_diagonally_dominant_problem(10000);
    std::cout<<"Solving..";
    auto x = jacobi(std::get<0>(tup), std::get<1>(tup), 1000);
    std::cout << are_ones(x);

}