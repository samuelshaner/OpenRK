/**
 * @file constants.h
 * @brief General library of constants.
 * @date May 10, 2015.
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

/** The states of the simulation */
#define NUM_STATES 8
#define START 0
#define PREVIOUS_OUT 1
#define PREVIOUS_IN 2
#define CURRENT 3
#define FORWARD_IN_OLD 4
#define FORWARD_OUT 5
#define FORWARD_OUT_OLD 6
#define END 7

/** The surfaces of a cube */
#define NUM_SURFACES 6
#define SURFACE_X_MIN 0
#define SURFACE_X_MAX 1
#define SURFACE_Y_MIN 2
#define SURFACE_Y_MAX 3
#define SURFACE_Z_MIN 4
#define SURFACE_Z_MAX 5

/** The boundary conditions */
#define VACUUM 0
#define REFLECTIVE 1

/** The transient methods */
#define IQS 0
#define ADIABATIC 1
#define THETA 2

/** The maximum number of iterations allowed for a power method eigenvalue
 *  solve in linalg.cpp */
#define MAX_LINALG_POWER_ITERATIONS 25000

/** The maximum number of iterations allowed for a linear solve in linalg.cpp */
#define MAX_LINEAR_SOLVE_ITERATIONS 1000

/** The tolerance on the initial and adiabatic flux solves */
#define FLUX_SOLVE_TOLERANCE 1.e-6

#endif /* CONSTANTS_H_ */
