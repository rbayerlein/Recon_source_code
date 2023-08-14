/*!
 * \file petsys_timer.h
 *
 * \brief Define a Timer class to display time information
 *
 * \author Jian Zhou
 * \date 2011-10-05
 * \since 0.1
 *
 * Copyright (c) 2011, Jian Zhou. All rights reserved.
 */
#ifndef PETSYS_TIMER_H
#define PETSYS_TIMER_H

#include <iostream>
#include <ctime>
#include <cstdio>
#include <sys/time.h>
#include <petsys_log.h>

namespace UCDPETSYS
{

class Timer
{
public:
    Timer();
    ~Timer();
public:
    void start();
    void stop();
    static void print(); // show local time

private:
    time_t m_t0;
    time_t m_t1;
    timeval m_hres_t0;
    timeval m_hres_t1;
};

}



#endif
