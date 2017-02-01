/**
 * @file constants.h
 * @brief General library of constants.
 * @date May 10, 2015.
 * @author Samuel Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

/** The maximum number of iterations allowed for a power method eigenvalue
 *  solve in linalg.cpp */
#define MAX_LINALG_POWER_ITERATIONS 25000

/** The maximum number of iterations allowed for a linear solve in linalg.cpp */
#define MAX_LINEAR_SOLVE_ITERATIONS 1000

/** The tolerance on the initial and adiabatic flux solves */
#define FLUX_SOLVE_TOLERANCE 1.e-6

#endif /* CONSTANTS_H_ */
