#include "../include/common.hpp"
#include "../include/timer.hpp"
#include "../include/solvers.h"

using namespace dp;

int main(int argc, char *argv[]) {
//     nw should be used for internal vector operations despite the jacobi method run in parallel
    int nw = std::stoi(argv[1]);
    int problem_size = std::stoi(argv[2]);


    std::tuple<Matrix, Vector> tup;
    {
        timer t("\ngeneration");
        tup = generate_diagonally_dominant_problem(problem_size, nw);
    }

    const int n_iter = 50;

    Vector x;
    {
        timer t("sequential");
        x = jacobi_seq(std::get<0>(tup), std::get<1>(tup), n_iter, nw);
    }
    std::cout << are_ones(x) << "  " << std::endl;

    {
        timer t("native");
        x = jacobi_native(std::get<0>(tup), std::get<1>(tup), n_iter, nw);
    }

    std::cout << are_ones(x) << "  " << std::endl;

//    for (auto &it: x) std::cout << it << std::endl;
    return 0;
}