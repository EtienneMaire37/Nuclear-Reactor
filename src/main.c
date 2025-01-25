#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "reactor.h"

int main()
{
    uint32_t steps_per_second = 10000;
    reactor_t* reactor = reactor_create(
        1.01, 3.2e-7, 1e-4, 0.0065, 1e7, 1e5, 300.0, 1e5, 600.0, 1. / steps_per_second
    );

    for (uint16_t i = 0; i < 100; i++)
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