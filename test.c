/**
*	math_aero
*	Copyright (C) 2021  Norbert Kiszka <norbert at linux dot pl>
*	This program is free software; you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation; either version 2 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License along
*	with this program; if not, write to the Free Software Foundation, Inc.,
*	51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <libgen.h>
#include <errno.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include "math_aero.h"

int main(void)
{
	double a, b, c, d;
	
	double pressure_fl350;

	a = 100 * c_knots_to_metersPerSecond;
	
	printf("100 knt to mps: %.9lf\n", a);
	
	a = a * c_metersPerSecond_to_knots;
	
	printf("Back to knots: %.9lf\n", a);
	
	printf("Dry air gas constant (ISA, ICAO, Sea level, 101325 Pa 15 °C): %.9lf\n", c_rgas_stp_iso);
	
	printf("Water constant: %.9lf\n", c_rgas_water);
	
	a = calc_dew_point(273.15 + 15, 0.44);
	
	printf("Dew point at 15 °C and 44%% humidity: %.9lf K\n", a);
	
	printf("Dew point at 15 °C and 44%% humidity: %.9lf °C\n", a - 273.15);
	
	a = calc_air_saturation(273.15 + 15);
	
	printf("Air water saturation at 15 °C (Arden Buck equation): %.9lf Pa\n", a);
	
	a = calc_air_saturation_tetens(273.15 + 15);
	
	printf("Air water saturation at 15 °C (Tetens/Murray equation): %.9lf Pa\n", a);
	
	a = calc_air_gas_constant(101325, 273.15 + 15, 0);
	
	printf("Air gas constant at 101325 Pa, 15 °C and 0%% humidity: %.9lf\n", a);
	
	a = calc_air_gas_constant(101325, 273.15 + 35, 0.6);
	
	printf("Air gas constant at 101325 Pa, 35 °C and 60%% humidity: %.9lf\n", a);
	
	a = calc_air_pressure(101325, 0, 9.80665, 0);
	
	printf("Air pressure at sea level (QNH 101325 Pa, 15 °C, 1g = 9.80665 m/s^2, 0.65 K/100m): %.9lf Pa\n", a);
	
	a = calc_air_pressure(102000, 1000 * 0.3048, 9.805, -1);
	
	printf("Air pressure at 1000 ft MSL (QNH 1020 hPa, 14 °C, 1g = 9.805 m/s^2, 0.65 K/100m): %.9lf Pa\n", a);
	
	// 1 feet = 1 m
	pressure_fl350 = calc_air_pressure(101325, 35000 * 0.3048, 9.80665, 0);
	
	printf("Air pressure at FL350 (QNH 101325 Pa, 15 °C, 1g = 9.80665 m/s^2, 0.65 K/100m): %.9lf Pa (value used in later tests)\n", pressure_fl350);
	
	a = calc_air_pressure_mark_drela(101325, 35000 * 0.3048);
	
	printf("Air pressure at FL350 with Mark Drela method (QNH 101325 Pa): %.9lf Pa (value not used in later tests)\n", a);
	
	a = calc_air_density(101325, 273.15 + 15, 0);
	
	printf("Air density at 101325 Pa, 15 °C and 0%% humidity: %.9lf kg/m^3\n", a);
	
	b = 273.15 + 15 - (35000 * 0.3048 * 0.65 * 0.01);
	
	a = calc_air_density(pressure_fl350, b, 0);
	
	printf("Air density at FL350, 15 °C, 0.65 K/100m (%.2lf °C) and 0%% humidity: %.9lf kg/m^3\n", b - 273.15, a);
	
	b = 273.15 + 35 - (35000 * 0.3048 * 0.60 * 0.01);
	
	a = calc_air_density(pressure_fl350, b, 0);
	
	printf("Air density at FL350, 35 °C, 0.60 K/100m (%.2lf °C) and 0%% humidity: %.9lf kg/m^3\n", b - 273.15, a);
	
	a = calc_air_density(pressure_fl350, b, 0.95);
	
	printf("Air density at FL350, 35 °C, 0.60 K/100m (%.2lf °C) and 95%% humidity: %.9lf kg/m^3\n", b - 273.15, a);
	
	a = calc_dynamic_pressure(1.225, 101325, 100 * c_knots_to_metersPerSecond);
	
	printf("Dynamic pressure at 100 knt, 101325 Pa, 15 °C: %.9lf Pa\n", a);
	
	a = calc_drag(1.225, 122 * c_knots_to_metersPerSecond, 0.032, 172 * 0.3048 * 0.3048);
	
	printf("Cessna 172 Drag at 122 knt (CD = 0.032), 1.225 kg/m^3 (101325 Pa, 15 °C): %.9lf N\n", a);
	
	a = calc_cd(1.225, 122 * c_knots_to_metersPerSecond, a, 172 * 0.3048 * 0.3048);
	
	printf("Calculate CD from previous data: %.9lf\n", a);
	
	a = b = calc_lift(1.225, 122 * c_knots_to_metersPerSecond, 0.24, 172 * 0.3048 * 0.3048);
	
	printf("Calculate lift force with CL = 0.24, Rho = 1.225 kg/m^3, 122 knt, 172 ft^2: %.9lf N\n", a);
	
	a = calc_cl(1.225, 122 * c_knots_to_metersPerSecond, a, 172 * 0.3048 * 0.3048);
	
	printf("Calculate CL from previous data: %.9lf\n", a);
	
	a = calc_speed_for_force(1.225, 0.24, b, 172 * 0.3048 * 0.3048);
	
	printf("Calculate speed (TAS) needed for %.9lf N with 1.225 kg/m^3, CL = 0.24, 172 ft^2: %.9lf knt\n", b, a * c_metersPerSecond_to_knots);
	
	c = calc_speed_of_sound(calc_air_density(pressure_fl350, 273.15 + 15 - (35000 * 0.3048 * 0.65 * 0.01), 0), pressure_fl350);
	
	printf("Calculate LSS (local speed of sound) at FL350 with QNH 101325 Pa, 15 °C, 0%% humidity: %.9lf m/s, %.9f knt\n", c, c * c_metersPerSecond_to_knots);
	
	a = calc_speed_of_sound(calc_air_density(101325, 273.15 + 40, 0.99), 101325);
	
	printf("Calculate LSS (local speed of sound) at sea level with 101325 Pa, 40 °C, 99%% humidity: %.9lf m/s, %.9f knt\n", a, a * c_metersPerSecond_to_knots);
	
	a = calc_mach_number(450 * c_knots_to_metersPerSecond, c);
	
	printf("Calculate Mach speed at FL350, 450 knt TAS: %.9lf\n", a);
	
	a = calc_pitching_moment_airfoil(-0.05, 10 * 0.3048 * 0.3048, 5 * 0.3048 + 4 * 0.0254, 100 * c_knots_to_metersPerSecond, 1.225, 101325);
	
	printf("Calculate pitching moment (torque) in Cessna 172 wing root at CM = -0.05, 1.225 kg/m^3, 101325 Pa, 100 knt, 10 ft^3: %.9lf Nm\n", a);
	
	a = calc_pitching_moment_aircraft(12, 5000, -400, 15, 0, 0, 0, -777, 10, 30);
	
	printf("Calculate pitching moment (torque) in imaginary aircraft (CG: 12 m, Wings lift: 5000 N, Wings PM: -400 Nm, Wings Cp = 15 m, Lift fuselage: 0 N, PM fuselage: 0 Nm, Lift HS: -777 N, PM HS: 10 Nm, HS Cp = 30 m): %.9lf Nm\n", a);
	
	a = calc_dynamic_viscosity(273.15 + 15);
	
	printf("Calculate air dynamic viscosity at 15 °C: %.9lf Pa/s\n", a);
	
	a = calc_dynamic_viscosity(273.15 + 40);
	
	printf("Calculate air dynamic viscosity at 40 °C: %.9lf Pa/s\n", a);
	
	a = calc_dynamic_viscosity(273.15 - 40);
	
	printf("Calculate air dynamic viscosity at -40 °C: %.9lf Pa/s\n", a);
	
	//a = calc_barometric_altitude(101325, pressure_fl350);
	
	a = calc_barometric_altitude(101325, pressure_fl350);
	
	printf("Calculate barometric altitude (pressure altitude) at FL350 with previous calculated pressure (%.9lf Pa) and QNH 101325 Pa: %.9lf ft\n", pressure_fl350, a / 0.3048);
	
	b = 273.15 + 15 - (0.0065 * 35000 * 0.3048);
	
	//a = calc_density_altitude(101325, pressure_fl350, b, 9.80665, 0);
	
	a = calc_density_altitude(101325, pressure_fl350, b, 9.80665, 0);
	
	printf("Calculate density altitude with previous pressure (%.9lf Pa), QNH 101325 Pa, %.2lf °C (ISA 15 °C): %.9lf ft\n", pressure_fl350, b - 273.15, a / 0.3048);
	
	b = 273.15 + 35 - (0.65 * 35000 * 0.3048 * 0.01);
	
	a = calc_density_altitude(101500, pressure_fl350, b, 9.80665, 0);
	
	printf("Calculate density altitude with previous pressure (%.9lf Pa), QNH 101500 Pa, %.2lf °C (ISA + 20 °C = 35 °C), 9.80665 m/s^2, 0%% humidity: %.9lf ft\n", pressure_fl350, b - 273.15, a / 0.3048);
	
	a = calc_density_altitude(101500, pressure_fl350, b, 9.80665, 0.3);
	
	printf("Calculate density altitude with previous pressure (%.9lf Pa), QNH 101500 Pa, %.2lf °C (ISA + 20 °C = 35 °C), 9.80665 m/s^2, 30%% humidity: %.9lf ft\n", pressure_fl350, b - 273.15, a / 0.3048);
	
	b = 273.15 + 15 - (0.65 * 35000 * 0.3048 * 0.01);
	
	a = calc_temperature_ram_rise(b, 2);
	
	printf("Calculate temperature ram rise at FL350 (%.2lf K, %.2lf °C) and speed Mach 2.0: %.9lf\n", b, b - 273.15, a);
	
	c = calc_sat_from_ram_rise(a, 2);
	
	printf("Calculate SAT (outside air temperature) from previous calculated ram rise (%.2lf) and speed Mach 2.0: %.2lf K, %.2lf °C\n", a, c, c - 273.15);
	
	d = b + a; // TAT temperature
	
	printf("Calculate TAT from previous data (SAT + RR): %.2lf K, %.2lf °C\n", d, d - 273.15);
	
	c = calc_sat_from_tat(d, 2);
	
	printf("Calculate SAT from previous calculated TAT (%.2lf K) and speed Mach 2.0: %.2lf K, %.2lf °C\n", d, c, c - 273.15);
	
	a = calc_Mach_from_CAS(300 * c_knots_to_metersPerSecond, pressure_fl350);
	
	printf("Calculate Mach speed from CAS 300 knots at FL350 %.9lf\n", a);
	
	b = calc_air_density(pressure_fl350, 288.15 - (0.0065 * 35000 * 0.3048), 0);
	
	d = calc_speed_of_sound(b, pressure_fl350); // LSS
	
	a = calc_TAS_from_Mach(a, d);
	
	printf("Calculate TAS speed from previous Mach speed at FL350: %.9lf knt\n", a * c_metersPerSecond_to_knots);
	
	a = calc_EAS_from_TAS(a, b);
	
	printf("Calculate EAS speed from previous TAS speed at FL350: %.9lf knt\n", a * c_metersPerSecond_to_knots);
	
	a = calc_TAS_from_EAS(a, b);
	
	printf("Calculate TAS speed from previous EAS speed at FL350: %.9lf knt\n", a * c_metersPerSecond_to_knots);
	
	a = calc_CAS_from_TAS(a, pressure_fl350, d);
	
	printf("Calculate CAS speed from previous TAS speed at FL350: %.9lf knt\n", a * c_metersPerSecond_to_knots);
	
	mass_point_t *mps;
	mps = calloc(3, sizeof(mass_point_t));
	
	mps->distance = 1;
	mps->mass = 12;
	
	mps[sizeof(mass_point_t)].distance = 2;
	mps[sizeof(mass_point_t)].mass = 16.5;
	
	mps[sizeof(mass_point_t)*2].distance = 3;
	mps[sizeof(mass_point_t)*2].mass = 44;
	
	printf("Mass points for later calculations (distance from nose in meters / mass in kg):\n");
	
	printf("1: %.2lf / %.2lf\n", mps->distance, mps->mass);
	
	printf("2: %.2lf / %.2lf\n", (mps + sizeof(mass_point_t))->distance, (mps + sizeof(mass_point_t))->mass);
	
	printf("3: %.2lf / %.2lf\n", (mps + 2 * sizeof(mass_point_t))->distance, (mps + 2 * sizeof(mass_point_t))->mass);
	
	d = calc_cg(mps, 3);
	
	printf("Center of gravity: %.9lf m\n", d);
	
	a = calc_total_mass(mps, 3);
	
	printf("Total mass (sum of all points): %.9lf kg\n", a);
	
	a = calc_moment_of_inertia(d, mps, 3);
	
	printf("Moment of inertia at previous calculated CG: %.9lf kg/m^2\n", a);
	
	a = calc_angular_acceleration(10, mps, 3);
	
	printf("Angular acceleration with torque 10 Nm: %.9lf rad/s^2 (%.9lf °/s^2)\n", a, a * (180.0 / M_PI));
	
	printf("\nThank You for flying with test airlines!\n");
	
	return 0;
}
