#include "../include/task_pool.h"
#include "../include/solvers.h"
#include <functional>

dp::Vector jacobi_task_pool(const Matrix A, const Vector b, const int max_iter, int nw,
                                  std::function<bool(Vector &)> stopping_criteria) {
    int n = A.size();
    Matrix s = zeros(n,n);
    Vector x(n, 0);
    Vector x_new(n, 0);

    task_pool tp(nw);

    std::vector<int> chunks = define_chunks(nw, n);

    auto f = [&A, &s, &x, n] (int start, int end, int task_id) {
        for (int i = 0; i < n; i++) {
            for (int j = start; j < end; j++) {
                if(i!=j) {
                    s[i][task_id] += A[i][j] * x[j];
                }
            }
        }
    };



    for(int k=0; k<max_iter; k++){
        for(int w=0; w<nw; w++) {
            auto fx = std::bind(f,chunks[w],chunks[w+1],w);
            tp.submit(fx);
        }
    }
    tp.submit(nullptr);

    tp.start();

    return x;

}