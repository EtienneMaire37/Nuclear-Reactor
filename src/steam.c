#pragma once

void get_enthalpies(double pressure, double temperature, double* h_water, double* h_steam) 
{
    *h_water = 0.0;
    *h_steam = 0.0;
    
    double min_pressure_diff = 1e9;
    int best_index = 0;
    for (int i = 0; i < STEAM_TABLE_SIZE; i++) 
    {
        double diff = fabs(steam_table[i].pressure - pressure);
        if (diff < min_pressure_diff) 
        {
            min_pressure_diff = diff;
            best_index = i;
        }
    }
    
    steam_table_entry_t entry = steam_table[best_index];
    if (temperature <= entry.temperature) 
    {
        *h_water = entry.h_water;
        *h_steam = entry.h_steam;
    } 
    else 
    {
        int next_index = best_index + 1;
        if (next_index < STEAM_TABLE_SIZE && steam_table[next_index].pressure == entry.pressure) 
        {
            steam_table_entry_t next_entry = steam_table[next_index];
            double frac = (temperature - entry.temperature) / (next_entry.temperature - entry.temperature);
            *h_water = entry.h_water + frac * (next_entry.h_water - entry.h_water);
            *h_steam = entry.h_steam + frac * (next_entry.h_steam - entry.h_steam);
        } 
        else 
        {
            *h_water = entry.h_water;
            *h_steam = entry.h_steam;
        }
    }
}

double get_water_density(double pressure, double temperature, double boiling_point) 
{
    double h_water, h_steam;
    get_enthalpies(pressure, temperature, &h_water, &h_steam);
    
    return 1000 * (1 - 0.0002 * (temperature - boiling_point));  // ~-0.2% density/°C
}

double calculate_void_fraction(double h_liquid, double h_actual, double h_vapor) 
{
    if (h_actual <= h_liquid) return 0.0;
    if (h_actual >= h_vapor) return 1.0;
    return (h_actual - h_liquid) / (h_vapor - h_liquid);  // Steam quality approximation
}