#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "reactor.h"

int main()
{
    uint32_t steps_per_second = 1000;
    reactor_t* reactor = reactor_create(
    1.005,        // k
    3.2e-7,       // κ
    1e-4,         // Λ
    0.0065,       // β
    1e7,          // C = 1e7 J/K 
    1e5,          // h = 1e5 W/K
    300.0,        // Coolant temp
    1e4,          // Initial n
    300.0,        // Initial T
    -1e-5,        // doppler coefficient
    1. / steps_per_second // dt
    );


    for (uint16_t i = 0; i < 1000; i++)
    {
        printf("s : %d\n", i);
        printf("    n = %f\n", floor(reactor->n));
        printf("    T = %f K\n", reactor->T);
        printf("    k = %f\n", reactor->k);
        printf("    P = %f W\n", reactor->κ * reactor->n);
        for (uint16_t j = 0; j < steps_per_second; j++)
            reactor_step(reactor);
    }

    reactor_free(reactor);
}