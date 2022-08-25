#ifndef PARALLEL_JACOBI_TASK_POOL_H
#define PARALLEL_JACOBI_TASK_POOL_H
#include <iostream>
#include <thread>
#include <vector>
#include <deque>
#include <mutex>
#include <functional>
#include <condition_variable>

class task_pool
{
private:
    int nw;
    std::deque<std::function<void(int,int,int)>> tasks;
    std::condition_variable cond;
    std::mutex ll;
    std::vector<std::thread> tids;
    bool stop_flag = false;
public:
    task_pool(int nw);
    ~task_pool();

    bool body(int i);
    void submit(std::function<void(int,int,int)> f);
    void start();
    void stop();

};


#endif //PARALLEL_JACOBI_TASK_POOL_H
