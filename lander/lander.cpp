// Mars lander simulator
// Version 1.11
// Mechanical simulation functions
// Gabor Csanyi and Andrew Gee, August 2019

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

#include "lander.h"

void autopilot (void)
  // Autopilot to adjust the engine throttle, parachute and attitude control
{
  // Defines the key local variables
  double m; // Number variable: mass
  const double K_h = 0.005, K_p = 0.75;
  m = UNLOADED_LANDER_MASS + fuel * FUEL_CAPACITY * FUEL_DENSITY; 
  double h = position.abs() - MARS_RADIUS;
  double delta = (GRAVITY * MARS_MASS * m / position.abs2()) / MAX_THRUST;
  double P_out = - K_p * (0.5 + K_h * h + velocity * position.norm());

  // Throttle conditions

  if (P_out <= -delta) {
    throttle = 0;
  } else if (P_out < 1 - delta) {
    throttle = delta + P_out;
  } else {
    throttle = 1;
  }
}

void numerical_dynamics (void)
  // This is the function that performs the numerical integration to update the
  // lander's pose. The time step is delta_t (global variable).
{
  // Defines the key local variables
  const bool use_verlet = true; // True if Verlet, false if Eurler
  vector3d F_t, F_g, F_d, acceleration; // 3D vector variables: thrust, gravity and drag, acceleration
  double rho, m; // Number variables: air density, mass
  rho = atmospheric_density(position);

  // Mass
  m = UNLOADED_LANDER_MASS + fuel * FUEL_CAPACITY * FUEL_DENSITY; 

  // Forces and acceleration
  F_g = - GRAVITY * MARS_MASS * m / position.abs2() * position.norm();
  F_t = thrust_wrt_world();
  F_d = - 0.5 * rho * DRAG_COEF_LANDER * M_PI * pow(LANDER_SIZE, 2)  * velocity.abs2() * velocity.norm();
  
  if (parachute_status == DEPLOYED) {
    F_d += - 0.5 * rho * DRAG_COEF_CHUTE *  5 * pow(2*LANDER_SIZE, 2)  * velocity.abs2() * velocity.norm();
  }

  acceleration = (F_g + F_t + F_d) / m;

  // Integrator
  static vector3d previous_position;
  if (use_verlet != true or simulation_time == 0) {
    // Euler integrator
    previous_position = position;
    position = position + velocity * delta_t;
    velocity = velocity + acceleration * delta_t;
    
  } else {
    // Verlet integrator
    vector3d future_position;
    future_position = 2 * position - previous_position + acceleration * pow(delta_t, 2);
    velocity = (future_position - previous_position) / (2 * delta_t);
    previous_position = position;
    position = future_position;
    
  }
  

  // Here we can apply an autopilot to adjust the thrust, parachute and altitude
  if (autopilot_enabled) autopilot();

  // Here we can apply 3-axis stabilization to ensure the base is always pointing downwards
  if (stabilized_attitude) attitude_stabilization();
}

void initialize_simulation (void)
  // Lander pose initialization - selects one of 10 possible scenarios
{
  // The parameters to set are:
  // position - in Cartesian planetary coordinate system (m)
  // velocity - in Cartesian planetary coordinate system (m/s)
  // orientation - in lander coordinate system (xyz Euler angles, degrees)
  // delta_t - the simulation time step
  // boolean state variables - parachute_status, stabilized_attitude, autopilot_enabled
  // scenario_description - a descriptive string for the help screen

  scenario_description[0] = "circular orbit";
  scenario_description[1] = "descent from 10km";
  scenario_description[2] = "elliptical orbit, thrust changes orbital plane";
  scenario_description[3] = "polar launch at escape velocity (but drag prevents escape)";
  scenario_description[4] = "elliptical orbit that clips the atmosphere and decays";
  scenario_description[5] = "descent from 200km";
  scenario_description[6] = "";
  scenario_description[7] = "";
  scenario_description[8] = "";
  scenario_description[9] = "";

  switch (scenario) {

  case 0:
    // a circular equatorial orbit
    position = vector3d(1.2*MARS_RADIUS, 0.0, 0.0);
    velocity = vector3d(0.0, -3247.087385863725, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 1:
    // a descent from rest at 10km altitude
    position = vector3d(0.0, -(MARS_RADIUS + 10000.0), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 2:
    // an elliptical polar orbit
    position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
    velocity = vector3d(3500.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 3:
    // polar surface launch at escape velocity (but drag prevents escape)
    position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE/2.0);
    velocity = vector3d(0.0, 0.0, 5027.0);
    orientation = vector3d(0.0, 0.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 4:
    // an elliptical orbit that clips the atmosphere each time round, losing energy
    position = vector3d(0.0, 0.0, MARS_RADIUS + 100000.0);
    velocity = vector3d(4000.0, 0.0, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 5:
    // a descent from rest at the edge of the exosphere
    position = vector3d(0.0, -(MARS_RADIUS + EXOSPHERE), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 6:
    break;

  case 7:
    break;

  case 8:
    break;

  case 9:
    break;

  }
}
