#pragma once

typedef struct reactor 
{
    double k;   // k factor
    double k_control_rods;
    double kappa;   // Power proportionality constant
    double lambda;   // Mean generation time
    double beta;   // Delayed neutron fraction
    double C;   // Heat capacity
    double h;   // Heat transfer coefficient
    double alpha_doppler; // Doppler coefficient
    double efficiency;  // Thermal to electric conversion efficiency
    double boiling_point; // Coolant boiling point

    double lambda_i[6];
    double beta_i[6];
    double C_i[6];

    double n;               // Neutron population
    double T;               // Temperature
    double T_initial;       // Initial temperature
    double T_coolant;       // Coolant temperature
    double target_n;        // Target neutron population
    double P_thermal;       // Power output
    double P_electric;      // Electricity output

    double pid_kp;     // Proportional gain
    double pid_ki;     // Integral gain
    double pid_kd;     // Derivative gain
    double pid_error_sum;
    double pid_last_error;

    double X;        // Xenon-135 concentration
    double I;        // Iodine-135 concentration

    double dt;              // Time step
} reactor_t;

#define LAMBDA_I 2.9e-5     // I-135 decay constant (s^-1)
#define LAMBDA_X 2.1e-5     // Xe-135 decay constant (s^-1)
#define SIGMA_X  3e-18      // Xe-135 absorption cross section (cm^2)
#define GAMMA_I  0.064      // I-135 yield per fission
#define GAMMA_X  0.002      // Xe-135 yield per fission

#define CR_MIN   0.95       // Minimum control rod position
#define CR_MAX   1.05       // Maximum control rod position

reactor_t* reactor_create(
    double k, double power_proportionality_constant, double mean_generation_time,
    double delayed_neutron_fraction, double heat_capacity, double heat_transfer_coefficient,
    double coolant_temperature, double n, double temperature, double doppler_coefficient, 
    double target_n, double efficiency, double coolant_boiling_point, double dt)
{
    reactor_t* reactor = (reactor_t*)malloc(sizeof(reactor_t));
    reactor->k = reactor->k_control_rods = k;
    reactor->kappa = power_proportionality_constant;
    reactor->lambda = mean_generation_time;
    reactor->beta = delayed_neutron_fraction;
    reactor->C = heat_capacity;
    reactor->h = heat_transfer_coefficient;
    reactor->T_coolant = coolant_temperature;
    reactor->alpha_doppler = doppler_coefficient;

    reactor->X = 0;
    reactor->I = 0;
    
    reactor->lambda_i[0] = 0.0124;
    reactor->lambda_i[1] = 0.0305;
    reactor->lambda_i[2] = 0.111;
    reactor->lambda_i[3] = 0.301;
    reactor->lambda_i[4] = 1.14;
    reactor->lambda_i[5] = 3.01;

    reactor->beta_i[0] = 0.000215;
    reactor->beta_i[1] = 0.001424;
    reactor->beta_i[2] = 0.001274;
    reactor->beta_i[3] = 0.002568;
    reactor->beta_i[4] = 0.000748;
    reactor->beta_i[5] = 0.000273;

    memset(reactor->C_i, 0, sizeof(reactor->C_i));

    reactor->pid_kp = 1e-4; 
    reactor->pid_ki = 1e-5;
    reactor->pid_kd = 1e-6;
    reactor->pid_error_sum = 0;
    reactor->pid_last_error = 0;
    reactor->target_n = target_n;
    reactor->P_electric = 0;
    reactor->efficiency = efficiency;
    reactor->boiling_point = coolant_boiling_point;

    reactor->n = n;
    reactor->T = temperature;
    reactor->T_initial = reactor->T;
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
        (reactor->k * (1 - reactor->beta) - 1) * reactor->n / reactor->lambda
    );

    for (int i = 0; i < 6; i++)
    {
        reactor->n += reactor->dt * (
            reactor->lambda_i[i] * reactor->C_i[i]
        );
    }

    for (int i = 0; i < 6; i++)
    {
        reactor->C_i[i] += reactor->dt * (
            reactor->k * reactor->beta_i[i] * reactor->n / reactor->lambda - reactor->lambda_i[i] * reactor->C_i[i]
        );
    }

    double P = reactor->kappa * reactor->n; // Power output
    double Q = reactor->h * (reactor->T - reactor->T_coolant); // Heat loss

    reactor->P_thermal = P;

    reactor->T += reactor->dt * (
        (P - Q) / reactor->C
    );

    reactor->P_electric = reactor->T < reactor->boiling_point ? 0 : 
    tanh((reactor->T - reactor->boiling_point) / 4) * reactor->efficiency * reactor->P_thermal;

    double phi = reactor->n / reactor->lambda; // Neutron flux

    reactor->I += reactor->dt * (GAMMA_I * phi - LAMBDA_I * reactor->I);
    reactor->X += reactor->dt * (
        LAMBDA_I * reactor->I + GAMMA_X * phi - LAMBDA_X * reactor->X - SIGMA_X * phi * reactor->X
    );

    double delta_k_xenon = -SIGMA_X * reactor->X / 2.41e24;  // Assume Σ_f = 2.41e24 cm^-1
    double delta_k_doppler = reactor->alpha_doppler * (reactor->T - reactor->T_initial);
    double delta_k_mod = -2e-4 * (reactor->T - reactor->T_initial);

    double error = reactor->target_n - reactor->n;
    reactor->pid_error_sum += error * reactor->dt;
    double derivative = (error - reactor->pid_last_error) / reactor->dt;
    reactor->pid_last_error = error;

    double delta_k = reactor->pid_kp * error + reactor->pid_ki * reactor->pid_error_sum + reactor->pid_kd * derivative;

    reactor->k_control_rods += clamp(delta_k, -0.0001, 0.0001);
    reactor->k_control_rods = clamp(reactor->k_control_rods, CR_MIN, CR_MAX);

    reactor->k = reactor->k_control_rods + delta_k_doppler + delta_k_mod + delta_k_xenon;
}