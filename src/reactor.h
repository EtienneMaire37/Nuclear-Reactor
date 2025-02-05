#pragma once

typedef struct reactor 
{
    real_t k;   // k factor
    real_t k_control_rods;
    real_t kappa;   // Power proportionality constant
    real_t lambda;   // Mean generation time
    real_t beta;   // Delayed neutron fraction
    real_t C;   // Heat capacity
    real_t h;   // Heat transfer coefficient
    real_t alpha_doppler; // Doppler coefficient
    real_t efficiency;  // Thermal to electric conversion efficiency
    real_t boiling_point; // Coolant boiling point

    real_t mod_density;       // Current moderator density (kg/m³)
    real_t void_fraction;     // Steam void fraction (0-1)
    real_t ref_mod_density;   // Reference moderator density (kg/m³)

    real_t lambda_i[6];
    real_t beta_i[6];
    real_t C_i[6];

    real_t n;               // Neutron population
    real_t T;               // Temperature          (K)
    real_t T_initial;       // Initial temperature  (K)
    real_t T_coolant;       // Coolant temperature  (K)
    real_t target_n;        // Target neutron population
    real_t P_thermal;       // Power output       (W)
    real_t P_electric;      // Electricity output (W)
    real_t pressure;        // Core pressure (MPa)

    real_t pid_kp;     // Proportional gain
    real_t pid_ki;     // Integral gain
    real_t pid_kd;     // Derivative gain
    real_t pid_error_sum;
    real_t pid_last_error;

    real_t X;        // Xenon-135 concentration
    real_t I;        // Iodine-135 concentration

    real_t dt;              // Time step
} reactor_t;

#define MIN_NEUTRONS 1e-10
#define MAX_XENON    1e25    // Prevent overflow

#define LAMBDA_I     2.9e-5  // I-135 decay (s⁻¹)
#define LAMBDA_X     2.1e-5  // Xe-135 decay (s⁻¹)
#define SIGMA_X      3e-18   // Xe-135 absorption cross-section (cm²)
#define GAMMA_I      0.064   // I-135 yield
#define GAMMA_X      0.001   // Xe-135 yield
#define SIGMA_F      0.0438  // Macroscopic fission cross-section (cm⁻¹)
#define SIGMA_I      0.1     // I-135 absorption cross-section (cm²)

#define CR_MIN   0.95       // Minimum control rod position
#define CR_MAX   1.05       // Maximum control rod position

reactor_t* reactor_create(
    real_t k, real_t power_proportionality_constant, real_t mean_generation_time,
    real_t delayed_neutron_fraction, real_t heat_capacity, real_t heat_transfer_coefficient,
    real_t coolant_temperature, real_t n, real_t temperature, real_t doppler_coefficient, 
    real_t target_n, real_t efficiency, real_t coolant_boiling_point, real_t pressure, real_t dt);
void reactor_free(reactor_t* reactor);
void reactor_step(reactor_t* reactor);