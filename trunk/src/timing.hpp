#ifndef timing_hpp
#define timing_hpp

#ifndef WIN32
#include "sys/time.h"
#include <iostream>
using std::cout;

class timing
{
private:
    struct timeval timer_;
    struct timezone timezone_;

public:
    timing()
    {}

    inline void tick(double& time_counter)
    {
        gettimeofday(&timer_, &timezone_);
        time_counter -= 1000 * timer_.tv_sec + timer_.tv_usec / 1000;
    }

    inline void tack(double& time_counter)
    {
        gettimeofday(&timer_, &timezone_);
        time_counter += 1000 * timer_.tv_sec + timer_.tv_usec / 1000;
    }

    void reportTimeNeeded(const char* message, double time_counter)
    {
        cout << message << " | " << time_counter << " milliseconds |\n" ;
    }
};

#endif

#ifdef WIN32

#include <time.h>
#include <iostream>
using std::cout;

class timing
{
private:
    double timePassed_;
    long begin_;
    long end_;

public:
    timing()
    {}

    inline void tick(double& time_counter)
    {
        time_counter -= clock();
    }

    inline void tack(double& time_counter)
    {
        time_counter += clock();
    }

    void reportTimeNeeded(const char* message, double time_counter)
    {
        cout << message << " | " << time_counter << " milliseconds |\n" ;
    }
};

#endif
#endif

