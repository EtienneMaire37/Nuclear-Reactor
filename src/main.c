#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <wchar.h>

// #include <io.h>
// #include <fcntl.h>
#ifndef _O_U16TEXT
  #define _O_U16TEXT 0x20000
#endif

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))
#define clamp(t, a, b) min(max(t, a), b)

#include "real.h"

#include "steam.h"
#include "reactor.h"

#include "steam.c"
#include "reactor.c"

#define TARGET_N 1e16

int main()
{
    uint32_t steps_per_second = 10000;
    reactor_t* reactor = reactor_create(    // GE BWR/4
    CR_MIN,         // k
    3.04e-7,        // κ
    1e-4,           // Λ
    0.0065,         // β
    1e7,            // C 
    4e7,            // h
    215. + 273.15,  // Coolant temp
    1e4,            // Initial n
    215. + 273.15,  // Initial T
    -1e-5,          // doppler coefficient
    1e5,            // target n
    784. / 3000.,   // efficiency
    285. + 273.15,  // coolant boiling point
    7,              // 7MPa
    1. / steps_per_second // dt
    );

    wprintf(L"time,neutron population,core temperature,multiplication factor k,thermal power,xenon concentration,iodine concentration,electric power,control rods' positions\n");

    for (uint32_t i = 0; i < 100000; i++)
    {
        wprintf(L"%d,%f,%f,%f,%f,%f,%f,%f,%f\n", 
        i, reactor->n, reactor->T - 273.15, reactor->k, 
        reactor->P_thermal, reactor->X, reactor->I, reactor->P_electric, 
        1 - (reactor->k_control_rods - CR_MIN) / (CR_MAX - CR_MIN));

        for (uint32_t j = 0; j < steps_per_second; j++)
            reactor_step(reactor);

        if (reactor->target_n < TARGET_N)
            reactor->target_n *= 1.0003;
        else
            reactor->target_n = TARGET_N;
    }

    reactor_free(reactor);
}