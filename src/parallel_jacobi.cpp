#include "../include/solvers.h"
#include <functional>
#include <barrier>
#include <numeric>


Vector jacobi_native::solve(Matrix A, Vector b, const int max_iter, int nw,
                            std::function<bool(Vector &)> stopping_criteria) {
    int n = A.size();
    // setup threads
    std::vector<std::thread> tids(nw);

    Vector diag(n);
    for (int i = 0; i < n; i++) {
        diag[i] = A[i][i];
        A[i][i] = 0;
    }

    Vector x(n, 0);
    Vector x_new(n);

    int k = 0;

    auto on_completion = [&k, &x, &x_new, &stopping_criteria, max_iter]() noexcept {
        // locking not needed here
        if (stopping_criteria != nullptr && stopping_criteria(x_new)) {
#ifdef PRINT_ITER
            std::cout << "in " << k << " iter" << std::endl;
#endif
            k = max_iter;
        }
        k++;
        x = x_new;
    };
    std::barrier sync_point(nw, on_completion);

    auto f = [&A, &b, &diag, n, max_iter, &k, &x, &x_new, &sync_point](int start, int end, int tid) {
        while (k < max_iter) {
            for (int i = start; i < end; i++) {
                double s = std::inner_product(A[i].begin(), A[i].end(), x.begin(), 0.0);
                x_new[i] = (b[i] - s) / diag[i];
            }

            sync_point.arrive_and_wait();
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

