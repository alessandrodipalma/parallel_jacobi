#include <iostream>
#include "../common/common.hpp"

using namespace dp;


Vector jacobi(Matrix A, Vector b, int max_iter, int nw=1) {

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
    Vector x = zeros(b.size());
//    Vector x_new(b.size());
    for (int k = 0; k < max_iter; k++) {
        for (int i = 0; i < A.size(); i++) {
//            Vector A1(i), x1(i), A2(A.size()-i), x2(A.size()-i);

//            copy(A[i].begin(), A[i].begin() + i, A1.begin());
//            copy(x.begin(), x.begin() + i, x1.begin());

//            copy(A[i].begin() + i + 1, A[i].end() , A2.begin());
//            copy(x.begin() + i + 1, x.end() , x2.begin());

//            double s1 = dot(A1, x1);
//            double s2 = dot(A2, x2);
            double s1 = dot(A[i], x, 0,i);
            double s2 = dot(A[i], x, i+1, x.size());

            x[i] = (b[i] - s1 - s2) / A[i][i];
        }
    }
    return x;

}


int main(int argc, char *argv[]) {
    // nw should be used for internal vector operations despite the jacobi method run in parallel
    int nw = std::stoi(argv[1]);
    int problem_size = std::stoi(argv[2]);

    std::tuple<Matrix, Vector> tup;
    {
        timer t("\ngeneration");
        tup = generate_diagonally_dominant_problem(problem_size, nw);
    }

    Vector x;
    {
        timer t("sequential");
        x = jacobi(std::get<0>(tup), std::get<1>(tup), 50, nw);
    }
    {
        timer t("product");
        x = jacobi_prod(std::get<0>(tup), std::get<1>(tup), 50, nw);
    }

    std::cout << are_ones(x) << "  " << x[0];

//    for (auto &it: x) std::cout << it << std::endl;
    return 0;
}