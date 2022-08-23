//
// Created by al on 19/08/2022.
//

#ifndef PARALLEL_JACOBI_SOLVERS_H
#define PARALLEL_JACOBI_SOLVERS_H

#include "common.hpp"

using namespace dp;

Vector jacobi_native(Matrix A, Vector b, int max_iter, int nw);

Vector jacobi_seq(Matrix A, Vector b, int max_iter, int nw);



#endif //PARALLEL_JACOBI_SOLVERS_H
