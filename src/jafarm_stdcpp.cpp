#include "../include/solvers.h"
#include <functional>
#include <barrier>
#include <numeric>
#include <sstream>
#include "../include/timer.hpp"

Vector jafarm_cpp::solve(Matrix A, Vector b, const int max_iter, int nw,
                         std::function<bool(Vector &)> stopping_criteria) {
    int n = A.size();



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
#ifdef MEASURE_ITERATES

        std::cout<<"Finished iterate k"<<std::endl;
        std::stringstream str;
        str << "jafarm_cpp time inside barrier";
        timer t(str.str());
#endif
        if (stopping_criteria(x_new)) {
#ifdef PRINT_IT
            std::cout << "in " << k << " iter" << std::endl;
#endif
            k = max_iter;
        } else {
            k++;
            x.swap(x_new);
        }
    };

    std::barrier sync_point(nw, on_completion);

    std::function f = [&A, &b, &diag, max_iter, &k, &x, &x_new, &sync_point](int start, int end, int tid) {
        while (k < max_iter) {
            {
#ifdef MEASURE_ITERATES
                std::stringstream str;
                str << "jafarm_cpp iterate, th" << tid;
                timer t(str.str());
#endif
                for (int i = start; i < end; i++) {

                    double s = std::inner_product(A[i].begin(), A[i].end(), x.begin(), 0.0);
                    x_new[i] = (b[i] - s) / diag[i];
                }
            }

            sync_point.arrive_and_wait();
        }
    };

    std::vector<std::thread> tids(nw);
    {

#ifdef MEASURE_ITERATES
        std::stringstream str;
        str << "jafarm_cpp setup cost" << std::endl;
        timer t(str.str());
#endif
        std::vector<int> chunks = define_chunks(nw, n);

        for (int i = 0; i < nw; i++) {
            tids[i] = std::thread(f, chunks[i], chunks[i + 1], i);
        }
    }
        for (int i = 0; i < nw; i++) {
            tids[i].join();
        }


    return x_new;
}

