#include "../include/solvers.h"
#include <functional>
#include <barrier>
#include <numeric>


Vector jacobi_alter_native::solve(Matrix A, Vector b, const int max_iter, int nw,
                            std::function<bool(Vector &)> stopping_criteria) {
    int n = A.size();

    std::vector<std::thread> tids(nw);

    Vector diag(n);
    for (int i = 0; i < n; i++) {
        diag[i] = A[i][i];
        A[i][i] = 0;
    }

    Vector x(n, 0);
    Vector x_new(n);
    Matrix sigma_parts = zeros(nw, n);

    int k = 0;

    auto on_completion = [&k, &x, &x_new, &stopping_criteria, &sigma_parts, max_iter, n, &b, &diag]() noexcept {
        for (int i=0;i<n; i++) {
            double s = std::accumulate(sigma_parts[i].begin(), sigma_parts[i].end(), 0.0);
            x_new[i] = (b[i] - s) / diag[i];
            std::fill(sigma_parts[i].begin(), sigma_parts[i].end(), 0.0);
        }

        if (stopping_criteria(x_new)) {
#ifdef PRINT_IT
            std::cout << "in " << k << " iter" << std::endl;
#endif
            k=max_iter;
        } else {
            k++;
            x.swap(x_new);
        }
    };

    std::barrier sync_point(nw, on_completion);

    std::function f = [&A, &b, &diag, &sigma_parts, max_iter, &k, &x, n, &x_new, &sync_point](int start, int end, int tid) {
        while (k < max_iter) {
            for (int i = start; i < end; i++) {
                double s = std::inner_product(A[i].begin(), A[i].end(), x.begin(), 0.0);
                x_new[i] = (b[i] - s) / diag[i];

                for(int j=0; j<n; j++) {
                    sigma_parts[j][tid] += x_new[i] * A[i][j];
                }
            }
            sync_point.arrive_and_wait();
        }
    };

    std::vector<int> chunks = define_chunks(nw, n);

    for (int i = 0; i < nw; i++) {
        tids[i] = std::thread(f, chunks[i], chunks[i+1], i);
    }

    for (int i = 0; i < nw; i++) {
        tids[i].join();
    }
    return x_new;
}

