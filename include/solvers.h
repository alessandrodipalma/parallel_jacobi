//
// Created by al on 19/08/2022.
//

#ifndef PARALLEL_JACOBI_SOLVERS_H
#define PARALLEL_JACOBI_SOLVERS_H

#include <functional>
#include <atomic>
#include "common.hpp"

using namespace dp;


inline std::vector<int> define_chunks(int nw, int n) {
    std::vector<int> chunks(nw + 1, 0);
    int chunk_size = n / nw;
    int chunk_reminder = n % nw;
    for (int i = 1; i < nw; i++) {
        chunks[i] = chunks[i - 1] + chunk_size + 1 * (chunk_reminder-- > 0);
    }
    chunks.back() = n;
    return chunks;

}


Vector jacobi_seq(Matrix A, Vector b, int max_iter, int nw);

Vector jacobi_seq_separate_iter
        (Matrix A, Vector b, int max_iter, std::function<bool(Vector &)> stopping_criteria = nullptr);

class solver {
public:
    virtual Vector solve(Matrix, Vector, int, int, std::function<bool(Vector &)>)=0;
    virtual std::string name() = 0;
    virtual ~solver() {};
};

class jafarm_cpp: public solver {
public:

    Vector solve(Matrix, Vector, int, int, std::function<bool(Vector &)>) override;
    std::string name() {return "c+";}

    ~jafarm_cpp() {};
};

class jafarm_ff: public solver {
public:
    Vector solve(Matrix, Vector, int, int, std::function<bool(Vector &)>) override;
    std::string name() {return "ff";}
    ~jafarm_ff() {};
};

class jafarm_omp: public solver {
public:
    Vector solve(Matrix, Vector, int, int, std::function<bool(Vector &)>) override;
    std::string name() {return "om";}
    ~jafarm_omp() {};
};

class jared_omp: public solver {
public:
    Vector solve(Matrix, Vector, int, int, std::function<bool(Vector &)>) override;
    std::string name() {return "al";}
    ~jared_omp() {};
};

class jared_cpp: public solver {
public:
    Vector solve(Matrix, Vector, int, int, std::function<bool(Vector &)>) override;
    std::string name() {return "al";}
    ~jared_cpp() {};
};

class jared_ff: public solver {
public:
    Vector solve(Matrix, Vector, int, int, std::function<bool(Vector &)>) override;
    std::string name() {return "af";}
    ~jared_ff() {};
};

class jacobi_2red: public solver {
public:
    Vector solve(Matrix, Vector, int, int, std::function<bool(Vector &)>) override;
    std::string name() {return "2r";}
    ~jacobi_2red() {};
};

#endif //PARALLEL_JACOBI_SOLVERS_H
