#include "../include/common.hpp"
#include "../include/timer.hpp"
#include "../include/solvers.h"

#include <iostream>
#include <sstream>

using namespace dp;

int main(int argc, char *argv[]) {
//     nw should be used for internal vector operations despite the jacobi method run in parallel
    unsigned const int nw = std::stoi(argv[1]);
    unsigned const int problem_size = std::stoi(argv[2]);
    // if 0, select the ones stopping criteria and set max_iter to 1e4, otherwise set a fixed number of
    unsigned int n_iter = std::stoi(argv[3]);

    std::tuple<Matrix, Vector> tup;
    {
        timer t("\ngeneration");
        tup = generate_diagonally_dominant_problem(problem_size, nw);
    }

    Vector x;
    {
        timer t("sequen");
        x = jacobi_seq_separate_iter(std::get<0>(tup), std::get<1>(tup), n_iter, nw);
    }
    std::cout << are_ones(x)  <<"  " << std::endl;
//    for (auto &it: x) std::cout << it << " ";

    std::stringstream s;
    s << "native with " << nw << " workers";

    std::function<bool(Vector&)> stopping;

    {
        timer t(s.str());
        x = jacobi_task_pool(std::get<0>(tup), std::get<1>(tup), n_iter, nw, are_ones);
    }
    std::cout << are_ones(x) << "  " << std::endl;
    for (auto &it: x) std::cout << it << " ";

    return 0;
}