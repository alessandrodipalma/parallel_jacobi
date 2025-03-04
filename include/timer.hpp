//
// Created by al on 19/08/2022.
//

#ifndef PARALLEL_JACOBI_TIMER_HPP
#define PARALLEL_JACOBI_TIMER_HPP

#include <iostream>
#include <chrono>


#define START(timename) auto timename = std::chrono::system_clock::now();
#define STOP(timename,elapsed)  auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - timename).count();


class timer {
    std::chrono::system_clock::time_point start;
    std::chrono::system_clock::time_point stop;
    std::string message;
    using usecs = std::chrono::microseconds;
    using msecs = std::chrono::milliseconds;

private:
    long * us_elapsed;

public:

    timer(const std::string m) : message(m),us_elapsed((long *)NULL) {
        start = std::chrono::system_clock::now();
    }

    timer(const std::string m, long * us) : message(m),us_elapsed(us) {
        start = std::chrono::system_clock::now();
    }

    ~timer() {
        stop =
                std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed =
                stop - start;
        auto musec =
                std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();

        std::cout << message << " computed in\t" << musec
                  << std::endl;
        if(us_elapsed != NULL)
            (*us_elapsed) = musec;
    }
};

#endif //PARALLEL_JACOBI_TIMER_HPP
