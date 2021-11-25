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


#ifndef MATH_AERO_H
#define MATH_AERO_H

extern const double c_knots_to_metersPerSecond;
extern const double c_metersPerSecond_to_knots;
extern const double c_rgas_stp_iso;
extern const double c_rgas_water;

inline double mps_to_knt(double mps)
{
	return mps * c_metersPerSecond_to_knots;
}

inline double knt_to_mps(double knt)
{
	return knt * c_knots_to_metersPerSecond;
}

typedef struct mass_point_s
{
	double distance; // distance in meters from nose tip (center of mass of one single object)
	double mass; // in kg
} mass_point_t;

/**
 * Calculate dew point with David Bolton method.
 * @param temp Temperature in K
 * @param humidity Humidity in number (0.5 = 50%)
 * @return Dew point in K
 */
double calc_dew_point(double temp, double humidity);

/**
 * Calculate air water saturation (Arden Buck equation).
 * @param temp Temperature in K
 * @return Saturation pressure in Pa
 */
double calc_air_saturation(double temp);

/**
 * Calculate air water saturation (Tetens/Murray equation).
 * @note Arden Buck equation is more precise (up to ~3%).
 * @param temp Temperature in K
 * @return Saturation pressure in Pa
 */
double calc_air_saturation_tetens(double temp);

/**
 * Calculate air gas constant in humid air
 * @param pressure Pressure in Pa
 * @param temp Temperature in K
 * @param humidity Humidity in number (0.5 = 50%)
 * @return m2/s2/K
 */
double calc_air_gas_constant(double pressure, double temp, double humidity);

/**
 * Calculate air density (Rho).
 * Humidity has bigger affect on density in higher temperatures - in those conditions air is less dense, so its important in calculations.
 * @param pressure Pressure in Pa
 * @param temp Temperature in K
 * @param humidity Humidity in number (0.5 = 50%)
 * @return Air density in kg/m^3
 */
double calc_air_density(double pressure, double temp, double humidity);

/**
 * Calculate air density (Rho). Method by Herman Wobus.
 * @param pressure Pressure in Pa
 * @param temp Temperature in K
 * @param humidity Humidity in number (0.5 = 50%)
 * @return Air density in kg/m^3
 */
double calc_air_density_herman_wobus(double pressure, double temp, double humidity);

/**
 * Calculate air dynamic pressure
 * @param air_density Air density in kg/m^3
 * @param air_pressure Air pressure in Pa
 * @param velocity Object velocity relative to air (m/s)
 * @return Dynamic pressure in Pa
 */
double calc_dynamic_pressure(double air_density, double air_pressure, double velocity);

/**
 * Calculate fluid (air) drag.
 * @param density Fluid density in kg/m^3
 * @param velocity Velocity in m/s
 * @param drag_coefficient Drag coefficient
 * @param area Area in m^2
 * @return Drag in N
 */
double calc_drag(double density, double velocity, double drag_coefficient, double area);

/**
 * Calculate drag coefficient.
 * @param density Fluid density in kg/m^3
 * @param velocity Velocity in m/s
 * @param drag_coefficient Drag coefficient
 * @param area Area in m^2
 * @return Cd
 */
double calc_cd(double density, double velocity, double drag_force, double area);

/**
 * Calculate fluid (air) lift.
 * @param density Fluid density in kg/m^3
 * @param velocity Velocity in m/s
 * @param lift_coefficient Lift coefficient
 * @param area Area in m^2
 * @return Lift in N
 */
double calc_lift(double density, double velocity, double lift_coefficient, double area);

/**
 * Calculate lift coefficient.
 * @param density Fluid density in kg/m^3
 * @param velocity Velocity in m/s
 * @param lift_coefficient Lift coefficient
 * @param area Area in m^2
 * @return Cd
 */
double calc_cl(double density, double velocity, double lift_force, double area);

/**
 * Calculate speed for given force, lift coefficient etc.
 * @param density Fluid density in kg/m^3
 * @param lift_coefficient Lift coefficient
 * @param force lift force (mass * 1g = mass * 9.80665)
 * @param area Area in m^2
 * @return Speed in m/s
 */
double calc_speed_for_force(double density, double lift_coefficient, double force, double area);

/**
 * Calculate speed of sound in air via air density and air pressure
 * @param air_density Air density in kg/m^3
 * @param air_pressure Air pressure in Pa
 * @return Speed of sound in m/s
 */
double calc_speed_of_sound(double air_density, double air_pressure);

/**
 * Does anybody need this very comlicated mathematic function?
 * @param velocity Velocity in m/s
 * @param speed_of_sound Speed of sound in m/s
 * @return Mach number (Mach speed)
 */
double calc_mach_number(double velocity, double speed_of_sound);

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
double calc_pitching_moment_airfoil(double moment_coefficient, double area, double chord, double velocity, double density, double pressure);

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
double calc_pitching_moment_aircraft(double cg, double lift_wings, double pm_wings, double cp_wings, double lift_fuselage, double pm_fuselage, double cp_fuselage, double lift_hs, double pm_hs, double cp_hs);

/**
 * Air dynamic viscosity from Sutherlands Equation.
 * @param temp Temperature in K
 * @return Dynamic viscosity in Pa/s
 */
double calc_dynamic_viscosity(double temp);

/**
 * Calculate Reynolds number
 * @param velocity Velocity in m/s
 * @param chord Chord length in m
 * @param density Fluid density
 * @param temp Fluid temperature in K
 * @return Reynolds number
 */
double calc_reynold(double velocity, double chord, double density, double temp);

/**
 * Calculate Reynolds number (unsigned long)
 * @param velocity Velocity in m/s
 * @param chord Chord length in m
 * @param density Fluid density
 * @param temp Fluid temperature in K
 * @return Reynolds number
 */
unsigned long calc_reynold_ul(double velocity, double chord, double density, double temp);

/**
 * Calculate air pressure at given altitude
 * @note This function is very accurate, but up to altitude 11 000 m (36 089 ft). Above this altitude (36 089 ft), this cannot be safely used in any aircraft, because temperature lapse rate is much different than in troposphere. In other words: above troposhpere, this function will give error, that will be increasing with every feet.
 * @param qnh QNH in Pa
 * @param altitude Altitude in m
 * @param gravity_acceleration Local gravity acceleration in m/s^2 (use 9.80665 for FL or if unsure)
 * @param temp_deviation Temperature deviation from ISA (Example: sea level is 20 °C, then use 5. Example2: sea level is 11 °C, then use -4)
 * @return Pressure in Pa
 */
double calc_air_pressure(double qnh, double altitude, double gravity_acceleration, double temp_deviation);

/**
 * Calculate air pressure at given altitude - method probably by Mark Drela.
 * @note This funtion was not fully tested. Please dont use this function in any aircraft.
 * @param qnh QNH in Pa
 * @param altitude Altitude in m
 * @return Pressure in Pa
 */
double calc_air_pressure_mark_drela(double qnh, double altitude);

/**
 * Calculate barometric altitude (aka pressure altitude)
 * @note This function is very accurate, but up to altitude 11 000 m (36 089 ft). Above this altitude (36 089 ft), this cannot be safely used in any aircraft, because (above 36 089 ft) temperature lapse rate is much different than in troposphere. In other words: above troposhpere, this function will give error, that will be increased with every feet.
 * @param qnh QNH in Pa (use 101325 for determine FL)
 * @param pressure Pressure in Pa
 * @return altitude in m
 */
double calc_barometric_altitude(double qnh, double pressure);

/**
 * Calculate density altitude (true altitude)
 * @param qnh QNH in Pa
 * @param pressure Pressure in Pa
 * @param temp outside air temperature
 * @param gravity_acceleration local gravity_acceleration
 * @param humidity
 * @return altitude in m
 */
double calc_density_altitude(double qnh, double pressure, double temp, double gravity_acceleration, double humidity);

/**
 * Calculate temperature rise by air compressibility
 * @param sat outside air temperature
 * @param mach_speed object Mach speed
 * @return temperature delta in K
 */
double calc_temperature_ram_rise(double sat, double mach_speed);

/**
 * Calculate outside temperature from ram rise
 * @param rr Temperature ram rise (difference between TAT and SAT)
 * @param mach_speed object Mach speed
 * @return SAT temperature in K
 */
double calc_sat_from_ram_rise(double rr, double mach_speed);

/**
 * Calculate SAT from TAT with given Mach speed
 * @param tat TAT temperature
 * @param mach_speed Mach number
 * @return SAT (outside air temperature)
 */
double calc_sat_from_tat(double tat, double mach_speed);

/**
 * Method borrowed (rewrited from JavaScript and modified) from Joachim K. Hochwarth method.
 * Source http://www.hochwarth.com/misc/AviationCalculator.html
 * Version: 1.8.2 (Joachim K. Hochwarth).
 * @param cas CAS speed
 * @param pressure Air static pressure
 * @return Mach number (Mach speed)
 */
double calc_Mach_from_CAS(double cas, double pressure);

/**
 * Calculate TAS from Mach
 * @param mach Mach speed
 * @param a Local speed of sound (around whole aircraft)
 * @return TAS speed
 */
double calc_TAS_from_Mach(double mach, double a);

/**
 * Calculate EAS from TAS
 * @param tas TAS speed
 * @param density Local air density
 * @return EAS speed
 */
double calc_EAS_from_TAS(double tas, double density);

/**
 * Calculate TAS from EAS
 * @param eas EAS speed
 * @param density Local air density
 * @return TAS speed
 */
double calc_TAS_from_EAS(double eas, double density);

/**
 * Method borrowed (rewrited and modified from JavaScript) from Joachim K. Hochwarth method.
 * Source http://www.hochwarth.com/misc/AviationCalculator.html
 * Version: 1.8.2 (Joachim K. Hochwarth).
 * @param cas CAS speed
 * @param pressure Air static pressure
 * @param a Local speed of sound (around whole aircraft)
 * @return CAS speed
 */
double calc_CAS_from_TAS(double tas, double pressure, double a);

/**
 * Calculate center of gravity
 * @note distance (struct mass_point_t) cant be zero or less than zero (distance = distance from aircraft nose tip).
 * @param points Points array
 * @param n Number of points
 * @return CG position
 */
double calc_cg(mass_point_t *points, size_t n);

/**
 * Calculate sum of mass in points array
 * @param points Points array
 * @param n Number of points
 * @return Total mass
 */
double calc_total_mass(mass_point_t *points, size_t n);

/**
 * Calculate moment of inertia from mass points relative to CG (2D object - two axis)
 * @param cg Center of Gravity
 * @param points Points array
 * @param n Number of points
 * @return Moment of inertia in kg/m^2
 */
double calc_moment_of_inertia(double cg, mass_point_t *points, size_t n);

/**
 * Calculate (theoretical) angular acceleration with given points array
 * @param torque Torque in Nm
 * @param points Points array
 * @param n Number of points
 * @return Angular acceleration in rad/s^2
 */
double calc_angular_acceleration(double torque, mass_point_t *points, size_t n);

#endif
