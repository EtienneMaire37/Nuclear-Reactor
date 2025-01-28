#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))
#define clamp(t, a, b) min(max(t, a), b)

#include "reactor.h"

#define TARGET_N 9.6e15

int main()
{
    uint32_t steps_per_second = 1000;
    reactor_t* reactor = reactor_create(    // BWR-4
    CR_MIN,       // k
    3.2e-7,       // κ
    1e-4,         // Λ
    0.0065,       // β
    1e7,          // C 
    4e7,          // h
    215. + 273.15,         // Coolant temp
    1e4,          // Initial n
    215. + 273.15,         // Initial T
    -1e-5,        // doppler coefficient
    1e5,          // target n
    1. / steps_per_second // dt
    );


    for (uint32_t i = 0; i < 100000; i++)
    {
        printf("s : %d\n", i);
        printf("    n = %f\n", floor(reactor->n));
        printf("    T = %f K\n", reactor->T);
        printf("    k = %f\n", reactor->k);
        printf("    P = %f W\n", reactor->kappa * reactor->n);
        printf("    control_rods = %f\n", 1 - (reactor->k_control_rods - CR_MIN) / (CR_MAX - CR_MIN));
        for (uint32_t j = 0; j < steps_per_second; j++)
            reactor_step(reactor);

        if (reactor->target_n < TARGET_N)
            reactor->target_n *= 1.003;
        else
            reactor->target_n = TARGET_N;
    }

    reactor_free(reactor);
}