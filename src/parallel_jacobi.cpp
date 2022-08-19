#include "../include/solvers.h"
#include <thread>

dp::Vector jacobi_native(const Matrix A, const Vector b, const int max_iter, int nw= 1) {
    int n = A.size();
    // setup threads
    std::vector<std::thread> tids(nw);

    for (int i=0; i<nw; i++) {
        tids[i] = std::thread();
    }

    Vector x(n, 0);

    // define task
    auto f = [A,b,n,&x] (int start, int end, int iter_id) {
        for (int i=start; i<end; i++) {
            double s = 0;
            for (int j=0; j<n; j++) {
                if(j!=i)
                    s = s + A[i][j] * x[j];
            }
            x[i] = (b[i] - s) / A[i][i];
        }
    };

    // define boundaries
    std::vector<int> chunks(nw+1,0);
    int chunk_size = n / nw ;
    int chunk_reminder = n % nw;
    for (int i=1;i<nw; i++){
        chunks[i] = chunks[i-1] + chunk_size + 1*(chunk_reminder-- > 0);
    }
    chunks.back() = n;

    // create a queue for iterations.
    for (int k=0; k<max_iter; k++) {
        for (int i = 0; i<nw; i++){
            f(chunks[i], chunks[i+1], k);
        }
        // send tasks to queue

    }

    return x;

}