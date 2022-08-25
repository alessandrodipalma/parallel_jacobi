//
// Created by al on 24/08/2022.
//

#include "../include/task_pool.h"


task_pool::task_pool(int nw)
{

    this->nw = nw;
    tids.reserve(nw);
    for(int i=0; i<nw; i++){
        tids[i] = std::thread(&task_pool::body, this, i);
    }

}

task_pool::~task_pool()
{
    std::cout<<"task_pool created"<<std::endl;
}

void task_pool::submit(std::function<void(int,int,int)> f) {
    {
        std::unique_lock<std::mutex> lock(ll);
        tasks.push_back(f);
        std::cout << tasks.size() <<" task pushed" << std::endl;
    }
    cond.notify_one();
}


bool task_pool::body(int i) {
    while(true) {
        std::function<void(int,int,int)> t= [](int, int, int) {};
        {
            std::unique_lock<std::mutex> lock(ll);
            cond.wait(lock, [&] () {return(!tasks.empty()) || stop_flag; });
            if (!tasks.empty()) {
                t = tasks.front();
                tasks.pop_front();
                if(t == nullptr) {
                    stop();
                }
            }
            if (stop_flag)
                return stop_flag;
        }
        t(0,0,0);
    }}


void task_pool::start(){
    for(int i=0; i<nw; i++){
        tids[i].join();
    }
}

void task_pool::stop() {
    {
        std::unique_lock<std::mutex> lock(ll);
        stop_flag = true;
    }
    cond.notify_all();
}


