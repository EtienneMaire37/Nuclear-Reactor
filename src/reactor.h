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

    double mod_density;       // Current moderator density (kg/m³)
    double void_fraction;     // Steam void fraction (0-1)
    double ref_mod_density;   // Reference moderator density (kg/m³)

    double lambda_i[6];
    double beta_i[6];
    double C_i[6];

    double n;               // Neutron population
    double T;               // Temperature          (K)
    double T_initial;       // Initial temperature  (K)
    double T_coolant;       // Coolant temperature  (K)
    double target_n;        // Target neutron population
    double P_thermal;       // Power output       (W)
    double P_electric;      // Electricity output (W)
    double pressure;        // Core pressure (MPa)

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
    double target_n, double efficiency, double coolant_boiling_point, double pressure, double dt);
void reactor_free(reactor_t* reactor);
void reactor_step(reactor_t* reactor);