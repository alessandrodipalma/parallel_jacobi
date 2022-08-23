#include "../include/solvers.h"
#include <thread>
#include <mutex>

std::vector<int> define_chunks(int nw, int n);

dp::Vector jacobi_native(const Matrix A, const Vector b, const int max_iter, int nw = 1) {
    int n = A.size();
    // setup threads
    std::vector<std::thread> tids(nw);


    Vector x(n, 0);
    Vector x_new(n, 0);

    // define task
    std::mutex ll;

    int k = 0;
    int w_count = 0;
    auto f = [A, b, n, max_iter, &k, &x, &ll, &w_count, nw](int start, int end, int iter_id) {
        while (k < max_iter) {
            for (int i = start; i < end; i++) {
                double s = 0;
                for (int j = 0; j < n; j++) {
                    if (j != i)
                        s = s + A[i][j] * x[j];
                }
                x[i] = (b[i] - s) / A[i][i];

            }
            {
                std::unique_lock<std::mutex> lock(ll);
                if (k>=max_iter) {
                    return;
                }
                if (are_ones(x)){
                    std::cout << "computed in " << k << " iterations"<<std::endl;
                    k=max_iter;
                    return;
                }

                if(++w_count == nw){
                    k++;
                    w_count=0;
                }
            }
        }
    };


    std::vector<int> chunks = define_chunks(nw, n);

    // create a queue for iterations.
    // send tasks to queue
    for (int i = 0; i < nw; i++) {
        tids[i] = std::thread(f, chunks[i], chunks[i + 1], k);
    }


    for (int i = 0; i < nw; i++) {
        tids[i].join();
    }

    return x;

}

std::vector<int> define_chunks(int nw, int n) {
    std::vector<int> chunks(nw + 1, 0);
    int chunk_size = n / nw;
    int chunk_reminder = n % nw;
    for (int i = 1; i < nw; i++) {
        chunks[i] = chunks[i - 1] + chunk_size + 1 * (chunk_reminder-- > 0);
    }
    chunks.back() = n;
    return chunks;
}
