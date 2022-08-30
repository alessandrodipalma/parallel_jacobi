#include "../include/common.hpp"
#include "../include/timer.hpp"
#include "../include/solvers.h"

#include <iostream>
#include <sstream>
#include <memory>
#include <list>
using namespace dp;

int main(int argc, char *argv[]) {
//     nw should be used for internal vector operations despite the jacobi method run in parallel
    unsigned const int min_w = std::stoi(argv[1]);
    unsigned const int max_w = std::stoi(argv[2]);
    unsigned const int problem_size = std::stoi(argv[3]);
    // if 0, select the ones stopping criteria and set max_iter to 1e4, otherwise set a fixed number of
    unsigned int n_iter = std::stoi(argv[4]);

    std::tuple<Matrix, Vector> tup;
    {
        std::stringstream s;
        s <<  "generation of " << problem_size << " sized problem";
        timer t(s.str());
        tup = generate_diagonally_dominant_problem(problem_size);
    }


    Vector x;
    {
        timer t("sq");
        x = jacobi_seq_separate_iter(std::get<0>(tup), std::get<1>(tup), n_iter, are_ones);
    }
    std::cout << are_ones(x) << "  " << std::endl;
//    for (auto &it: x) std::cout << it << " ";

    std::vector<std::unique_ptr<solver>> solvers;
    solvers.emplace_back(new jacobi_native{});
    solvers.emplace_back(new jacobi_ff{});
    solvers.emplace_back(new jacobi_omp{});

    for (int i = 0; i<solvers.size(); i++){
        for (unsigned int nw = min_w; nw <= max_w; nw++) {

            {
                std::stringstream s;
                s << solvers[i]->name() << " with " << nw << " workers";
                timer t(s.str());
                x = solvers[i]->solve(std::get<0>(tup), std::get<1>(tup), n_iter, nw, are_ones);
            }
#ifdef CHECK_ONES
            std::cout << are_ones(x) << "  " << std::endl;
#endif
        }
    }
}