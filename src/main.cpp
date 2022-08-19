#include "../include/solvers.hpp"

using namespace dp;

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
        x = jacobi_seq(std::get<0>(tup), std::get<1>(tup), 50, nw);
    }
    {
        timer t("product");
        x = jacobi_native(std::get<0>(tup), std::get<1>(tup), 50, nw);
    }

    std::cout << are_ones(x) << "  " << x[0];

//    for (auto &it: x) std::cout << it << std::endl;
    return 0;
}