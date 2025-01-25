#pragma once

typedef struct reactor 
{
    double k;   // k factor
    double k_control_rods;
    double κ;   // Power proportionality constant
    double Λ;   // Mean generation time
    double β;   // Delayed neutron fraction
    double C;   // Heat capacity
    double h;   // Heat transfer coefficient

    double λi[6];
    double βi[6];
    double Ci[6];

    double n;               // Neutron population
    double T;               // Temperature
    double T_initial;       // Initial temperature
    double T_coolant;       // Coolant temperature

    double dt;              // Time step
} reactor_t;

reactor_t* reactor_create(
    double k, double power_proportionality_constant, double mean_generation_time,
    double delayed_neutron_fraction, double heat_capacity, double heat_transfer_coefficient,
    double coolant_temperature, double n, double temperature, double dt)
{
    reactor_t* reactor = (reactor_t*)malloc(sizeof(reactor_t));
    reactor->k = reactor->k_control_rods = k;
    reactor->κ = power_proportionality_constant;
    reactor->Λ = mean_generation_time;
    reactor->β = delayed_neutron_fraction;
    reactor->C = heat_capacity;
    reactor->h = heat_transfer_coefficient;
    reactor->T_coolant = coolant_temperature;
    
    reactor->λi[0] = 0.0124;
    reactor->λi[1] = 0.0305;
    reactor->λi[2] = 0.111;
    reactor->λi[3] = 0.301;
    reactor->λi[4] = 1.14;
    reactor->λi[5] = 3.01;

    reactor->βi[0] = 0.000215;
    reactor->βi[1] = 0.001424;
    reactor->βi[2] = 0.001274;
    reactor->βi[3] = 0.002568;
    reactor->βi[4] = 0.000748;
    reactor->βi[5] = 0.000273;

    memset(reactor->Ci, 0, sizeof(reactor->Ci));

    reactor->n = n;
    reactor->T = temperature;
    reactor->T_initial = temperature;
    reactor->dt = dt;
    return reactor;
}

void reactor_free(reactor_t* reactor)
{
    free(reactor);
}

void reactor_step(reactor_t* reactor)
{
    reactor->n += reactor->dt * (
        (reactor->k * (1 - reactor->β) - 1) * reactor->n / reactor->Λ
    );

    for (int i = 0; i < 6; i++)
    {
        reactor->n += reactor->dt * (
            reactor->λi[i] * reactor->Ci[i]
        );
    }

    for (int i = 0; i < 6; i++)
    {
        reactor->Ci[i] += reactor->dt * (
            reactor->k * reactor->βi[i] * reactor->n / reactor->Λ - reactor->λi[i] * reactor->Ci[i]
        );
    }

    double P = reactor->κ * reactor->n; // Power output
    double Q = reactor->h * (reactor->T - reactor->T_coolant); // Heat loss

    reactor->T += reactor->dt * (
        (P - Q) / reactor->C
    );

    double delta_k_doppler = -1e-5 * (reactor->T - reactor->T_initial);
    reactor->k = reactor->k_control_rods + delta_k_doppler;
}