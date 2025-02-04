#pragma once

reactor_t* reactor_create(
    real_t k, real_t power_proportionality_constant, real_t mean_generation_time,
    real_t delayed_neutron_fraction, real_t heat_capacity, real_t heat_transfer_coefficient,
    real_t coolant_temperature, real_t n, real_t temperature, real_t doppler_coefficient, 
    real_t target_n, real_t efficiency, real_t coolant_boiling_point, real_t pressure, real_t dt)
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
    reactor->pressure = pressure;

    reactor->X = 1e13;
    reactor->I = 1e14;
    
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
    reactor->ref_mod_density = get_water_density(pressure, temperature - 273.15, coolant_boiling_point - 273.15);

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
    // real_t theta = 0.5;  // Implicit-explicit balance
    // real_t rho = reactor->k - 1.0;
    // real_t prompt_term = (rho * (1 - reactor->beta) - 1) / reactor->lambda;

    reactor->n += reactor->dt * (
        ((1 - reactor->beta) * reactor->k - 1) * reactor->n / reactor->lambda
    );
    reactor->n = clamp(reactor->n, MIN_NEUTRONS, 1e20);
    for (uint8_t i = 0; i < 6; i++) 
    {
        reactor->n += reactor->dt * reactor->lambda_i[i] * reactor->C_i[i];
        reactor->C_i[i] += reactor->dt * (
            reactor->beta_i[i] * reactor->k * reactor->n / reactor->lambda - reactor->lambda_i[i] * reactor->C_i[i]
        );
    }

    real_t P = reactor->kappa * reactor->n; // Power output
    real_t Q = reactor->h * (reactor->T - reactor->T_coolant); // Heat loss

    reactor->P_thermal = P;

    reactor->T += reactor->dt * (
        (P - Q) / reactor->C
    );

    real_t h_water, h_steam;
    get_enthalpies(reactor->pressure, reactor->T - 273.15, &h_water, &h_steam);  // T in °C

    reactor->mod_density = get_water_density(reactor->pressure, reactor->T - 273.15, reactor->boiling_point - 273.15);
    reactor->void_fraction = calculate_void_fraction(h_water, 
                                    h_water + (reactor->T - reactor->T_coolant) * 4.1868, 
                                    h_steam);

    h_water *= 1000;
    h_steam *= 1000;

    real_t density_ratio = reactor->ref_mod_density / reactor->mod_density;
    real_t delta_k_density = 0.003 * (density_ratio - 1);
    real_t delta_k_void = -0.015 * reactor->void_fraction;

    real_t delta_k_mod = delta_k_density + delta_k_void;

    real_t delta_h = h_steam - h_water;
    real_t m_dot = (delta_h > 0) ? reactor->P_thermal / delta_h : 0.;
    
    reactor->P_electric = m_dot * (h_steam - h_water) * reactor->efficiency * 0.95;

    if (reactor->P_electric < 0) reactor->P_electric = 0.;

    // real_t phi = reactor->n / reactor->lambda; // Neutron flux
    real_t v = 2.2e5;  // Thermal neutron velocity (cm/s)
    real_t phi = reactor->n / v / reactor->lambda;  // n/cm²/s

    // Implicit Integration
    real_t termI = GAMMA_I * SIGMA_F * phi;
    reactor->I = (reactor->I + reactor->dt * termI) / (1 + reactor->dt * LAMBDA_I);
    reactor->I = clamp(reactor->I, 0.0, 1e20);

    real_t termX = LAMBDA_I * reactor->I + GAMMA_X * SIGMA_F * phi;
    real_t denomX = 1 + reactor->dt * (LAMBDA_X + SIGMA_X * phi);
    reactor->X = (reactor->X + reactor->dt * termX) / denomX;
    reactor->X = clamp(reactor->X, 0.0, MAX_XENON);

    real_t delta_k_xenon = -SIGMA_X * reactor->X / SIGMA_F;
    real_t delta_k_doppler = reactor->alpha_doppler * (reactor->T - reactor->T_initial);

    real_t error = reactor->target_n - reactor->n;
    reactor->pid_error_sum += error * reactor->dt;
    reactor->pid_error_sum = clamp(reactor->pid_error_sum, -1e5, 1e5);  // Reduced integral windup
    
    real_t derivative = (error - reactor->pid_last_error) / reactor->dt;
    reactor->pid_last_error = error;

    real_t delta_k = reactor->pid_kp * error + 
                     reactor->pid_ki * reactor->pid_error_sum + 
                     reactor->pid_kd * derivative;
    
    reactor->k_control_rods += clamp(delta_k * reactor->dt, -0.001, 0.001);  // Smaller adjustments
    reactor->k_control_rods = clamp(reactor->k_control_rods, CR_MIN, CR_MAX);

    reactor->k = reactor->k_control_rods + delta_k_doppler + delta_k_mod + delta_k_xenon;
    // printf("%f %f %f\n", delta_k_doppler, delta_k_mod, delta_k_xenon);
}