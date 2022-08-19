#ifndef PARALLEL_JACOBI_SOLVERS_HPP
#define PARALLEL_JACOBI_SOLVERS_HPP

#include "common.hpp"

using namespace dp;

Vector jacobi_native(Matrix A, Vector b, int max_iter, int nw);
Vector jacobi_seq(Matrix A, Vector b, int max_iter, int nw);


#endif //PARALLEL_JACOBI_SOLVERS_HPP
