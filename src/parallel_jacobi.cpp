#include "../include/solvers.h"
#include <thread>
#include <mutex>
#include <functional>

dp::Vector jacobi_native(const Matrix A, const Vector b, const int max_iter, int nw,
                         std::function<bool(Vector &)> stopping_criteria) {
    int n = A.size();
    // setup threads
    std::vector<std::thread> tids(nw);


    Vector x(n, 0);

    // define task
    std::mutex ll;

    std::atomic<unsigned int> k = 0;
    std::atomic<unsigned int> w_count = 0;
    auto f = [&A, &b, n, max_iter, &k, &x, &w_count, nw, stopping_criteria](int start, int end, int tid) {
        while (k < max_iter) {
            for (int i = start; i < end; i++) {
                double s = 0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        s = s + A[i][j] * x[j];
                    }
                }
                x[i] = (b[i] - s) / A[i][i];

            }

            if (k >= max_iter) {
                return;
            }
            if (stopping_criteria != nullptr && stopping_criteria(x)) {
                std::cout << "computed in " << k << " iterations" << std::endl;
                k = max_iter; // to let the other workers return
                return;
            }
            if (++w_count == nw) {
                k++;
                w_count = 0;
            }
        }
    };


    std::vector<int> chunks = define_chunks(nw, n);

    // create a queue for iterations.
    // send tasks to queue
    for (int i = 0; i < nw; i++) {
        tids[i] = std::thread(f, chunks[i], chunks[i + 1], i);
    }


    for (int i = 0; i < nw; i++) {
        tids[i].join();
    }
    return x;
}

