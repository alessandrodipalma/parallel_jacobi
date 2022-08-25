//#include "../include/solvers.h"
//#include <thread>
//#include <mutex>
//#include <functional>
//#include <condition_variable>
//#include <deque>
//#include <jmorecfg.h>
//
//bool vector_and(std::vector<bool> x) {
//    bool res = true;
//    for (bool b: x) {
//        res = res && b;
//    }
//    return res;
//}
//
//
//dp::Vector jacobi_native_sep_iter(const Matrix A, const Vector b, const int max_iter, int nw,
//                                  std::function<bool(Vector &)> stopping_criteria) {
//    int n = A.size();
//    std::vector<int> chunks = define_chunks(nw, n);
//    std::vector<std::thread> tids(nw);
//
//    Vector x(n, 0);
//    Vector x_new(n, 0);
//    int k = 0;
//    int global_count = 0;
//
//    std::mutex mtx;
//    std::condition_variable cv;
//    std::vector<std::vector<std::deque<int>>> jqueue(nw,std::vector<std::deque<int>>(1));
//    for (int i = 0; i < nw; i++) {
//        for (int j = 0; j < nw; j++)
//            jqueue[i][0].push_back(j);
//    }
//
//    auto notify_available = [&mtx, &cv, &global_count, &k, nw, &jqueue, &x, &x_new](int chunk, int iter) {
//        std::unique_lock<std::mutex> lck(mtx);
//        for (int i = 0; i < nw; i++) {
//            jqueue[i][iter+1].push_back(chunk);
//        }
//        global_count++;
//        cv.notify_all();
//        std::cout << "computed chunk " << chunk << ", global_count" << global_count << std::endl;
//        if (global_count == nw) {
//            k++;
//            x = x_new;
//            global_count = 0;
////            std::cout<<"concluded iter "<<k<<std::endl;
//        };
//    };
//
//    auto f = [&A, &b, &cv, &mtx, nw, max_iter,
//            &k, &x, &x_new,
//            &chunks, &jqueue, &notify_available](int start, int end, int tid) {
//
//        while (k < max_iter) {
//            Vector s(end - start, 0);
//            std::vector<bool> computed_chunks(nw, false);
//             // store sigmas for every x[i] to be computed
//
//            while (!vector_and(computed_chunks)) {
//                std::unique_lock<std::mutex> lck(mtx);
//                while (jqueue[tid].empty()) {
//                    cv.wait(lck);
//                }
//                while (!jqueue[tid].empty()) {
//                    int c = jqueue[tid][k].front();
//                    jqueue[tid][k].pop_front();
//                    for (int i = start; i < end; i++) {
//                        for (int j = chunks[c]; j < chunks[c + 1]; j++) {
//                            if(i!=j) {
//                                s[i] += A[i][j] * x[j];
//                            }
//                        }
//                    }
//                    computed_chunks[c] = true;
//                }
//            }
////            std::cout<< tid<<"computed all chunks" <<std::endl;
//            for (int i = start; i < end; i++)
//                x_new[i] = (b[i] - s[i]) / A[i][i];
//
//            notify_available(tid, k);
//        }
//    };
//
//    // create a queue for iterations.
//    // send tasks to queue
//    for (int i = 0; i < nw; i++) {
//        tids[i] = std::thread(f, chunks[i], chunks[i + 1], i);
//    }
//
//    for (int i = 0; i < nw; i++) {
//        tids[i].join();
//    }
//
//    return x;
//}
//
