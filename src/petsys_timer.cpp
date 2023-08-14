/*
 *  petsys_timer.cpp
 *  ucd_pet_recon
 *
 *  Created by Jian Zhou on 12/7/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <petsys_timer.h>
using namespace UCDPETSYS;

Timer::Timer()
{
}

Timer::~Timer()
{
}

void Timer::start()
{
    time(&m_t0);
    gettimeofday(&m_hres_t0, NULL);
    SystemLog::write("Timer is started on %s", ctime(&m_t0));
}

void Timer::stop()
{
    time(&m_t1);
    gettimeofday(&m_hres_t1, NULL);
    SystemLog::write("Timer is stop on %s", ctime(&m_t1));

    double et = (m_hres_t1.tv_sec - m_hres_t0.tv_sec) * 1000.0;
    et += (m_hres_t1.tv_usec - m_hres_t0.tv_usec) / 1000.0;
    SystemLog::write("Elasped time is %.2lf ms.\n", et);
}

void Timer::print()
{
    time_t t;
    time(&t);
    SystemLog::write("%s", ctime(&t));
}
