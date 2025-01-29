#pragma once

typedef struct steam_table_entry
{
    double pressure;    // MPa
    double temperature; // °C
    double h_water;     // Specific enthalpy of liquid water (kJ/kg)
    double h_steam;     // Specific enthalpy of saturated steam (kJ/kg)
} steam_table_entry_t;

steam_table_entry_t steam_table[] = 
{
    // Pressure (MPa) | Temperature (°C) | h_water (kJ/kg) | h_steam (kJ/kg)
    {7.0, 250, 1087.3, 2799.5},
    {7.0, 275, 1185.2, 2772.8},
    {7.0, 285, 1230.4, 2756.1},  // Typical BWR-4 operating point
    {7.0, 300, 1315.8, 2725.3},
    {7.5, 285, 1220.7, 2748.9},
    {6.5, 285, 1240.1, 2763.2},
};

#define STEAM_TABLE_SIZE (sizeof(steam_table) / sizeof(steam_table_entry_t))