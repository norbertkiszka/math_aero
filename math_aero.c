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

// In many cases, multipling takes less CPU time than dividing.
const double c_knots_to_metersPerSecond  = 0.51444444444444444444444444444444444;
const double c_metersPerSecond_to_knots  = 1.94384449244060464323524684004951268;

const double c_rgas_stp_iso                   = 287.05287424704397380992304533720016479; // [J/kg K], ISA, ICAO (101325 Pa, 15 °C, 1.225 kg/1m^3)
const double c_rgas_water                     = 461.5; // Same as above, but for water

/**
 * Calculate dew point with David Bolton method.
 * @param temp Temperature in K
 * @param humidity Humidity in number (0.5 = 50%)
 * @return Dew point in K
 */
double calc_dew_point(double temp, double humidity)
{
	if(humidity <= 0)
		return 0; // :)
	
	if(humidity > 1) // wet wet wet
		humidity = 1;
	
	temp -= 273.15;
	double gamma_T_RH = log(humidity * exp((17.67 - (temp / 243.5)) * (temp / (243.5 + temp))));
	return ((243.5 * gamma_T_RH) / (17.67 - gamma_T_RH))+273.15;
}

/**
 * Calculate air water saturation (Arden Buck equation).
 * @param temp Temperature in K
 * @return Saturation pressure in Pa
 */
double calc_air_saturation(double temp)
{
	if(temp < 273.15)
		return 611.15 * exp((23.036-(temp-273.15)/333.7)*((temp-273.15)/(temp+6.67)));

	return 611.21 * exp((18.678-(temp-273.15)/234.5)*((temp-273.15)/(temp-16.01)));
}

/**
 * Calculate air water saturation (Tetens/Murray equation).
 * @note Arden Buck equation is more precise (up to ~3%).
 * @param temp Temperature in K
 * @return Saturation pressure in Pa
 */
double calc_air_saturation_tetens(double temp)
{
	if(temp < 273.15)
		return 610.78 * exp((21.875*(temp-273.15))/(temp-7.65));

	return 610.78 * exp((17.27*(temp-273.15))/(temp-35.85));
}

/**
 * Calculate air gas constant in humid air
 * @param pressure Pressure in Pa
 * @param temp Temperature in K
 * @param humidity Humidity in number (0.5 = 50%)
 * @return m2/s2/K
 */
double calc_air_gas_constant(double pressure, double temp, double humidity)
{
	// safety checks
	if(humidity < 0)
		humidity = 0;
	if(humidity > 1)
		humidity = 1;
	
	//double psat = calc_air_saturation(temp);
	//return c_rgas_stp_iso / (1 - (humidity * psat/pressure) * (1 - c_rgas_stp_iso/c_rgas_water));
	
	//return c_rgas_stp_iso / (1 - (humidity * calc_air_saturation(temp)/pressure) * (1 - c_rgas_stp_iso/c_rgas_water));
	
	return c_rgas_stp_iso / (1 - (humidity * calc_air_saturation(temp)/pressure) * 0.37800027248744538788116642535896972);
	
	
}

/**
 * Calculate air density (Rho).
 * Humidity has bigger affect on density in higher temperatures - in those conditions air is less dense, so its important in calculations.
 * @param pressure Pressure in Pa
 * @param temp Temperature in K
 * @param humidity Humidity in number (0.5 = 50%)
 * @return Air density in kg/m^3
 */
double calc_air_density(double pressure, double temp, double humidity)
{
	// safety checks
	if(humidity < 0)
		humidity = 0;
	if(humidity > 1)
		humidity = 1;
	
	double rh = humidity > 0 ? calc_air_gas_constant(pressure, temp, humidity) : c_rgas_stp_iso;
	
	// Calculate and return value:
	return pressure / (rh * temp);
}

/**
 * Calculate air density (Rho). Method by Herman Wobus.
 * @param pressure Pressure in Pa
 * @param temp Temperature in K
 * @param humidity Humidity in number (0.5 = 50%)
 * @return Air density in kg/m^3
 */
double calc_air_density_herman_wobus(double pressure, double temp, double humidity)
{
	// safety checks
	if(humidity < 0)
		humidity = 0;
	if(humidity > 1)
		humidity = 1;
	
	double dew = calc_dew_point(temp, humidity) - 273.15;
	double p = 0.99999683 + dew*(-0.90826951E-02 + dew*(0.78736169E-04 + dew*(-0.61117958E-06 + dew*(0.43884187E-08 + dew*(-0.29883885E-10 + dew*(0.21874425E-12 + dew*(-0.17892321E-14 + dew*(0.11112018E-16 + dew*(-0.30994571E-19)))))))));
	

    // Step 2: calculate the vapor pressure.
	
	//double psat_mbar = 6.1078 / (pow(p, 8));
	//double pv_pascals = psat_mbar * 100.0;
	double pv_pascals = 610.78 / (pow(p, 8));

    // Step 3: calculate the pressure of dry air, given the vapor
    // pressure and actual pressure.
	
	double pd_pascals = pressure - pv_pascals;
	
    // Step 4: calculate the air density, using the equation for
    // the density of a mixture of dry air and water vapor.
	
	return
        (pd_pascals / (c_rgas_stp_iso * temp)) +
        (pv_pascals / (c_rgas_water * temp));
}

/**
 * Calculate air dynamic pressure
 * @param air_density Air density in kg/m^3
 * @param air_pressure Air pressure in Pa
 * @param velocity Object velocity relative to air (m/s)
 * @return Dynamic pressure in Pa
 */
double calc_dynamic_pressure(double air_density, double air_pressure, double velocity)
{
	//double mach_number = calc_mach_number(velocity, calc_speed_of_sound(air_density, air_pressure));
	//return mach_number * mach_number * 0.5 * 1.4 * air_pressure;

	return 0.5 * air_density * velocity * velocity;
}

/**
 * Calculate fluid (air) drag.
 * @param density Fluid density in kg/m^3
 * @param velocity Velocity in m/s
 * @param drag_coefficient Drag coefficient
 * @param area Area in m^2
 * @return Drag in N
 */
double calc_drag(double density, double velocity, double drag_coefficient, double area)
{
	return 0.5 * drag_coefficient * density * velocity * velocity * area;
}

/**
 * Calculate drag coefficient.
 * @param density Fluid density in kg/m^3
 * @param velocity Velocity in m/s
 * @param drag_coefficient Drag coefficient
 * @param area Area in m^2
 * @return Cd
 */
double calc_cd(double density, double velocity, double drag_force, double area)
{
	return 2 * drag_force / (area * density * velocity * velocity);
}

/**
 * Calculate fluid (air) lift.
 * @param density Fluid density in kg/m^3
 * @param velocity Velocity in m/s
 * @param lift_coefficient Lift coefficient
 * @param area Area in m^2
 * @return Lift in N
 */
double calc_lift(double density, double velocity, double lift_coefficient, double area)
{
	return lift_coefficient * area * density * velocity * velocity * 0.5;
}

/**
 * Calculate lift coefficient.
 * @param density Fluid density in kg/m^3
 * @param velocity Velocity in m/s
 * @param lift_coefficient Lift coefficient
 * @param area Area in m^2
 * @return Cd
 */
double calc_cl(double density, double velocity, double lift_force, double area)
{
	return 2 * lift_force / (area * density * velocity * velocity);
}

/**
 * Calculate speed for given force, lift coefficient etc.
 * @param density Fluid density in kg/m^3
 * @param lift_coefficient Lift coefficient
 * @param force lift force (mass * 1g = mass * 9.80665)
 * @param area Area in m^2
 * @return Speed in m/s
 */
double calc_speed_for_force(double density, double lift_coefficient, double force, double area)
{
	return sqrt(force / (0.5 * area * density * lift_coefficient));
}

/**
 * Calculate speed of sound in air via air density and air pressure
 * @param air_density Air density in kg/m^3
 * @param air_pressure Air pressure in Pa
 * @return Speed of sound in m/s
 */
double calc_speed_of_sound(double air_density, double air_pressure)
{
	return sqrt((air_pressure*1.4)/air_density);
}

/**
 * Does anybody need this very comlicated mathematic function?
 * @param velocity Velocity in m/s
 * @param speed_of_sound Speed of sound in m/s
 * @return Mach number (Mach speed)
 */
double calc_mach_number(double velocity, double speed_of_sound)
{
	return velocity / speed_of_sound;
}

/**
 * Calculate moment (torque) generated by airfoil
 * @param moment_coefficient Airfoil moment coefficient
 * @param area Airfoil area in m^2
 * @param chord Chord length in m
 * @param velocity Velocity in m/s
 * @param density Fluid density in kg/m^3
 * @param pressure Fluid pressure in Pa
 * @return Torque in Nm
 */
double calc_pitching_moment_airfoil(double moment_coefficient, double area, double chord, double velocity, double density, double pressure)
{
	//return moment_coefficient * calc_dynamic_pressure(air_density, air_pressure, velocity) * area * chord;
	return moment_coefficient * 0.5 * density * velocity * velocity * area * chord;
}

/**
 * Calculate pitching moment (torque) of whole aircraft in relative to cg
 * @param cg center of gravity
 * @param lift_wings Lift force (25% chord) generated by both wings
 * @param pm_wings Moment (torque at 25% chord) generated by wings
 * @param cp_wings Center of pressure in wings (25% chord in most cases)
 * @param lift_fuselage Lift force (25% chord) generated by fuselage
 * @param pm_fuselage Moment (torque at 25% chord) generated by fuselage
 * @param cp_fuselage Center of pressure in fuselage (25% chord in most cases)
 * @param lift_hs Lift force (25% chord) generated by horizontal stabilizer(s)
 * @param pm_hs Moment (torque at 25% chord) generated by horizontal stabilizer
 * @param cp_wings Center of pressure in horizontal stabilizer(s) (25% chord in most cases)
 * @return Torque in Nm
 */
double calc_pitching_moment_aircraft(double cg, double lift_wings, double pm_wings, double cp_wings, double lift_fuselage, double pm_fuselage, double cp_fuselage, double lift_hs, double pm_hs, double cp_hs)
{
	return
		lift_wings * (cg - cp_wings) + pm_wings * (cg - cp_wings) +
		lift_fuselage * (cg - cp_fuselage) + pm_fuselage * (cg - cp_fuselage) + 
		lift_hs * (cg - cp_hs) + pm_hs * (cg - cp_hs);
}

/**
 * Air dynamic viscosity from Sutherlands Equation.
 * @param temp Temperature in K
 * @return Dynamic viscosity in Pa/s
 */
double calc_dynamic_viscosity(double temp)
{
	// Sea level air viscosity (15 °C): 0.0000181206
	return 0.0000181206 * pow(temp / 288.15, 1.5) * (288.15+113) / (temp+113);
}

/**
 * Calculate Reynolds number
 * @param velocity Velocity in m/s
 * @param chord Chord length in m
 * @param density Fluid density
 * @param temp Fluid temperature in K
 * @return Reynolds number
 */
double calc_reynold(double velocity, double chord, double density, double temp)
{
	return density * velocity * chord / calc_dynamic_viscosity(temp);
}

/**
 * Calculate Reynolds number (unsigned long)
 * @param velocity Velocity in m/s
 * @param chord Chord length in m
 * @param density Fluid density
 * @param temp Fluid temperature in K
 * @return Reynolds number
 */
unsigned long calc_reynold_ul(double velocity, double chord, double density, double temp)
{
	return (unsigned long)(density * velocity * chord / calc_dynamic_viscosity(temp));
}

/**
 * Calculate air pressure at given altitude
 * @note This function is very accurate, but up to altitude 11 000 m (36 089 ft). Above this altitude (36 089 ft), this cannot be safely used in any aircraft, because temperature lapse rate is much different than in troposphere. In other words: above troposhpere, this function will give error, that will be increasing with every feet.
 * @param qnh QNH in Pa
 * @param altitude Altitude in m
 * @param gravity_acceleration Local gravity acceleration in m/s^2 (use 9.80665 for FL or if unsure)
 * @param temp_deviation Temperature deviation from ISA (Example: sea level is 20 °C, then use 5. Example2: sea level is 11 °C, then use -4)
 * @return Pressure in Pa
 */
double calc_air_pressure(double qnh, double altitude, double gravity_acceleration, double temp_deviation)
{
	//return qnh * pow(((288.15-0.0065*altitude)/288.15), (gravity_acceleration*0.0289644)/(8.3144598*0.0065));
	//return qnh * pow(((288.15-0.0065*altitude)/288.15), gravity_acceleration/(c_rgas_stp_iso*0.0065));
	//return qnh * pow(((288.15-0.0065*altitude)/288.15), gravity_acceleration/1.86584368260578559173268331505823880);
	
	return qnh * pow((((288.15+temp_deviation)-0.0065*altitude)/(288.15+temp_deviation)), gravity_acceleration/1.86584368260578559173268331505823880);
}

/**
 * Calculate air pressure at given altitude - method probably by Mark Drela.
 * @note This funtion was not fully tested. Please dont use this function in any aircraft.
 * @param qnh QNH in Pa
 * @param altitude Altitude in m
 * @return Pressure in Pa
 */
double calc_air_pressure_mark_drela(double qnh, double altitude)
{
	const double altkm = altitude * 0.001;
	return qnh*exp(-0.118*altkm-(0.0015*altkm*altkm)/(1-0.018*altkm+0.0011*altkm*altkm));
}

/**
 * Calculate barometric altitude (aka pressure altitude)
 * @note This function is very accurate, but up to altitude 11 000 m (36 089 ft). Above this altitude (36 089 ft), this cannot be safely used in any aircraft, because (above 36 089 ft) temperature lapse rate is much different than in troposphere. In other words: above troposhpere, this function will give error, that will be increased with every feet.
 * @param qnh QNH in Pa (use 101325 for determine FL)
 * @param pressure Pressure in Pa
 * @return altitude in m
 */
double calc_barometric_altitude(double qnh, double pressure)
{
	//return (pow(pressure / qnh, (8.3144598*0.0065)/(gravity_acceleration*0.0289644)) * 288.15 - 288.15) / -0.0065;
	//return (pow(pressure / qnh, (c_rgas_stp_iso*0.0065)/(gravity_acceleration)) * 288.15 - 288.15) / -0.0065;
	//return (pow(pressure / qnh, (1.86584368260578559173268331505823880)/(gravity_acceleration)) * 288.15 - 288.15) / -0.0065;
	return (pow(pressure / qnh, 0.19026310540355634293163689108041581) * 288.15 - 288.15) / -0.0065;
}

/**
 * Calculate density altitude (true altitude)
 * @param qnh QNH in Pa
 * @param pressure Pressure in Pa
 * @param temp outside air temperature
 * @param gravity_acceleration local gravity_acceleration
 * @param humidity
 * @return altitude in m
 */
double calc_density_altitude(double qnh, double pressure, double temp, double gravity_acceleration, double humidity)
{
	double sea_level_temperature = temp + calc_barometric_altitude(101325, pressure) * 0.0065;
	
	return (pow(pressure / qnh, (calc_air_gas_constant(pressure, temp, humidity)*0.0065)/gravity_acceleration) * sea_level_temperature - sea_level_temperature) / -0.0065;
}

/**
 * Calculate temperature rise by air compressibility
 * @param sat outside air temperature
 * @param mach_speed object Mach speed
 * @return temperature delta in K
 */
double calc_temperature_ram_rise(double sat, double mach_speed)
{
	return sat * 0.2 * mach_speed * mach_speed;
}

/**
 * Calculate outside temperature from ram rise
 * @param rr Temperature ram rise (difference between TAT and SAT)
 * @param mach_speed object Mach speed
 * @return SAT temperature in K
 */
double calc_sat_from_ram_rise(double rr, double mach_speed)
{
	return rr / (0.2 * mach_speed * mach_speed);
}

/**
 * Calculate SAT from TAT with given Mach speed
 * @param tat TAT temperature
 * @param mach_speed Mach number
 * @return SAT (outside air temperature)
 */
double calc_sat_from_tat(double tat, double mach_speed)
{
	return tat / (0.2 * mach_speed * mach_speed + 1);
}

/**
 * Method borrowed (rewrited from JavaScript and modified) from Joachim K. Hochwarth method.
 * Source http://www.hochwarth.com/misc/AviationCalculator.html
 * Version: 1.8.2 (Joachim K. Hochwarth).
 * @param cas CAS speed
 * @param pressure Air static pressure
 * @return Mach number (Mach speed)
 */
double calc_Mach_from_CAS(double cas, double pressure)
{
	//double CaSLSI = sqrt(CGamma*CRGasSI*288.15);
	//double CaSLSI = 340.29399054347109199181176109050284140;
	
	//return 2.23606797749978969641442005933384962*sqrt(pow(((101325/pressure)*(pow((cas*cas/(5.0*CaSLSI*CaSLSI)) + 1, (CGamma/(CGamma - 1))) - 1) + 1), (CGamma - 1)/CGamma) - 1);
	//return 2.23606797749978969641442005933384962*sqrt(pow(((101325/pressure)*(pow((cas*cas/(5.0*CaSLSI*CaSLSI)) + 1, (1.4/(1.4 - 1))) - 1) + 1), (1.4 - 1)/1.4) - 1);
	//return 2.23606797749978969641442005933384962*sqrt(pow(((101325/pressure)*(pow((cas*cas/578999.99999999988358467817306518554687500) + 1, (1.4/(1.4 - 1))) - 1) + 1), (1.4 - 1)/1.4) - 1);
	return 2.23606797749978969641442005933384962*sqrt(pow(((101325/pressure)*(pow((cas*cas/578999.99999999988358467817306518554687500) + 1, 3.5) - 1) + 1), 0.28571428571428564291423413123993669) - 1);
}

/**
 * Calculate TAS from Mach
 * @param mach Mach speed
 * @param a Local speed of sound (around whole aircraft)
 * @return TAS speed
 */
double calc_TAS_from_Mach(double mach, double a)
{
	return mach * a;
}

/**
 * Calculate EAS from TAS
 * @param tas TAS speed
 * @param density Local air density
 * @return EAS speed
 */
double calc_EAS_from_TAS(double tas, double density)
{
	return tas * sqrt(density/1.225);
}

/**
 * Calculate TAS from EAS
 * @param eas EAS speed
 * @param density Local air density
 * @return TAS speed
 */
double calc_TAS_from_EAS(double eas, double density)
{
	return eas * sqrt(1.225/density);
}

/**
 * Method borrowed (rewrited and modified from JavaScript) from Joachim K. Hochwarth method.
 * Source http://www.hochwarth.com/misc/AviationCalculator.html
 * Version: 1.8.2 (Joachim K. Hochwarth).
 * @param cas CAS speed
 * @param pressure Air static pressure
 * @param a Local speed of sound (around whole aircraft)
 * @return CAS speed
 */
double calc_CAS_from_TAS(double tas, double pressure, double a)
{
	//return sqrt(5)*CaSLSI*sqrt(pow(((pressure/101325)*(pow((tas*tas/(5.0*a*a)) + 1, (CGamma/(CGamma - 1))) - 1) + 1), (CGamma - 1)/CGamma) - 1);
	//return 760.92049518987191208951870180499099661*sqrt(pow(((pressure/101325)*(pow((tas*tas/(5.0*a*a)) + 1, (CGamma/(CGamma - 1))) - 1) + 1), (CGamma - 1)/CGamma) - 1);
	//return 760.92049518987191208951870180499099661*sqrt(pow(((pressure/101325)*(pow((tas*tas/(5.0*a*a)) + 1, (1.4/(1.4 - 1))) - 1) + 1), (1.4 - 1)/1.4) - 1);
	return 760.92049518987191208951870180499099661*sqrt(pow(((pressure/101325)*(pow((tas*tas/(5.0*a*a)) + 1, 3.5) - 1) + 1), 0.28571428571428564291423413123993669) - 1);
}

/**
 * Calculate center of gravity
 * @note distance (struct mass_point_t) cant be zero or less than zero (distance = distance from aircraft nose tip).
 * @param points Points array
 * @param n Number of points
 * @return CG position
 */
double calc_cg(mass_point_t *points, size_t n)
{
	size_t i;
	mass_point_t *p;
	double dividend = 0;
	double divisor = 0;
	
	for(i = 0; i < n; i++)
	{
		p = points + i * sizeof(mass_point_t);
		dividend += p->distance * p->mass;
		divisor += p->mass;
	}
	
	return dividend / divisor;
}

/**
 * Calculate sum of mass in points array
 * @param points Points array
 * @param n Number of points
 * @return Total mass
 */
double calc_total_mass(mass_point_t *points, size_t n)
{
	size_t i;
	double mass = 0;
	
	for(i = 0; i < n; i++)
	{
		mass += (points + i * sizeof(mass_point_t))->mass;
	}
	
	return mass;
}

/**
 * Calculate moment of inertia from mass points relative to CG (2D object - two axis)
 * @param cg Center of Gravity
 * @param points Points array
 * @param n Number of points
 * @return Moment of inertia in kg/m^2
 */
double calc_moment_of_inertia(double cg, mass_point_t *points, size_t n)
{
	size_t i;
	double dist;
	double moment_of_inertia = 0;
	
	for(i = 0; i < n; i++)
	{
		dist = fabs((points + i * sizeof(mass_point_t))->distance - cg);
		moment_of_inertia += (points + i * sizeof(mass_point_t))->mass * dist * dist;
	}
	
	return moment_of_inertia;
}

/**
 * Calculate (theoretical) angular acceleration with given points array
 * @param torque Torque in Nm
 * @param points Points array
 * @param n Number of points
 * @return Angular acceleration in rad/s^2
 */
double calc_angular_acceleration(double torque, mass_point_t *points, size_t n)
{
	double cg = calc_cg(points, n);
	double moment_of_inertia = calc_moment_of_inertia(cg, points, n);
	
	return torque / moment_of_inertia;
}
