///
/// @file
///
/// @brief Contains all the global constants, as pre-processor definitions, used in the code.
///
/// This MPCD code avoids the use of explicit numbers where ever possible, and instead makes use of pre-processor defines. These are
/// stated in this file, with comments explaining each.
/// Broadly speaking, these are broken up into the following categories:
/// - Compile options
/// - Program constants
/// - Mathematical universal constants
/// - Dimensionality (_1D, _2D, _3D)
/// - Mersenne Twister random number generator constants
/// - Global flags (output, streaming, particle initial conditions, nematic options, hydrodynamic interactions, non-ideal EOS, multiphase)
/// - Thermostat options
/// - Collision operators
/// - Boundary condition types
/// - MD coupling options
/// - Monte Carlo options
/// - Swimmer types
/// - Debug (log verbosity) options
///
/// Note that as flags are commented in/out, the Doxygen output is not reliable for this file. See comments in the code
/// for more information.
///

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ GLOBAL PARAMETERS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

/* ****************************************** */
/* ************ COMPILE OPTIONS ************* */
/* ****************************************** */
/// @brief Enable debugging/ logging to terminal. Comment out to disable. Verbosity level is controlled by the
/// `"debugOut"` parameter in the json file.
# define DBG
/// @brief Force flush the output buffer, enforcing files to be printed immediately. Comment out to disable.
# define FFLSH
// /// @brief Use the legacy Mersenne Twister RNG. Comment out to use the Xoroshiro RNG.
// # define RNG_MERSENNE
/// @brief Enable the floating point signal handlers.
// #define FPE

/* ****************************************** */
/* ************ PROGRAM CONSTANTS *********** */
/* ****************************************** */
/// @brief Number of bins used for distributions. Best if an odd integer.
# define BINS 101
/// @brief Tolerance for the precision of the code. Used in: invert2x2(), checkplacement(),checkplace(),chooseT(),secant_time().
# define TOL 1e-8
/// @brief Tolerance for active rotations in activeSRD().
# define ROTTOL 1E-5
/// @brief Small times have a magnitude less than SCALE.
# define SCALE 1E-2
/// @brief Maximum number of times a particle can bounce between boundary conditions in a single time step.
# define NBOUNCE 50
/// @brief Maximum number of allowed MPCD species.
# define MAXSPECI 10
/// @brief Maximum number of allowed MPCD boundary conditions.
# define MAXBC 50
/// @brief Maximum number of vertices allowed for boundary conditions in bcMode==VERT_SURF
# define MAXVERT 10

/* ****************************************** */
/* *********** UNIVERSAL CONSTANTS ********** */
/* ****************************************** */
/// @brief Pi.
# define pi M_PI
/// @brief Euler's number.
# define e M_E

/* ****************************************** */
/* ************* DIMENSIONALITY ************* */
/* ****************************************** */
/// @brief Constant for 3D.
# define _3D 3
/// @brief Constant for 2D.
# define _2D 2
/// @brief Constant for 1D.
# define _1D 1

/* ****************************************** */
/* ********* MERSENNE RANDOM NUMBER ********* */
/* ****************************************** */
/// @brief Word count for Mersenne Twister random number generation.
# define NN 624
/// @brief Middle word for Mersenne Twister random number generation.
# define MM 397
/// @brief Constant vector a for Mersenne Twister random number generation.
# define MATRIX_A 0x9908b0dfUL
/// @brief Most significant w-r bits for Mersenne Twister random number generation.
# define UPPER_MASK 0x80000000UL
/// @brief Least significant r bits for Mersenne Twister random number generation.
# define LOWER_MASK 0x7fffffffUL

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* *************** GLOBAL FLAGS ************* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

/* ****************************************** */
/* ***************** OUTPUT ***************** */
/* ****************************************** */
/// @brief If set of data should be printed out to a file set flag to OUT.
# define OUT 1
/// @brief If set of data should NOT be printed, set flag to NOUT.
# define NOUT 0

/* ****************************************** */
/* *************** FILE NAMES *************** */
/* ****************************************** */
/// @brief Assumed string length for file names, etc. Long because clusters often have long paths
# define STRLN 500

/* ****************************************** */
/* ***************** STREAM ***************** */
/* ****************************************** */
/// @brief If the particle needs to stream its STREAM.
# define STREAM 1
/// @brief If the particle does not need to stream its NO_STREAM.
# define NO_STREAM 0

/* ****************************************** */
/* ******************* QDIST **************** */
/* ****************************************** */
//The following global variables are flag values for the variable QDIST in species structures.
//Used for both MPCD particles and swimmers --- qDist in json.
/// @brief Particle position distribution. PPF indicates that the particles can be placed randomly and freely within the control volume.
# define PRF 0
/// @brief Particle position distribution. READ indicates that the particles' positions should be read from a file.
#define READ 1
/// @brief Particle position distribution. CENTRED places a swimmer at the centre of the control volume.
#define CENTRED 2

/* ****************************************** */
/* ******************* VDIST **************** */
/* ****************************************** */
//The following global variables are flag values for the variable VDIST in species structures.
//Used for both MPCD particles and swimmers --- vDist in json.
/// @brief Velocity distribution. RANDVEL indicates that the particles' velocities are drawn from a have uniformly random velocity distribution.
# define RANDVEL 0
// /// @brief Velocity distribution. READ must be applied to both position and velocity.
// #define READ 1
/// @brief Velocity distribution. HOMO indicates that the particles all have the same initial speed but different directions.
# define HOMO 2
/// @brief Velocity distribution. Indicates that the particles' initial velocity is given by a given value along each axis in either direction.
# define HOMOAXIS 3
/// @brief Velocity distribution. GAUSS indicates that the particles' velocity follow a gaussian distribution.
# define GAUSS 4

/* ****************************************** */
/* ******************* UDIST **************** */
/* ****************************************** */
//The following global variables are flag values for the variable UDIST in species structures.
//Used for both MPCD particles and swimmers --- oDist in json.
/// @brief Orientation distribution. RANDORIENT indicates that the particles' direction are drawn from a have uniformly random direction.
# define RANDORIENT 0
/// @brief Orientation distribution. Align all particles along X axis.
# define ALIGNX 1
/// @brief Orientation distribution. Align all particles along Y axis.
# define ALIGNY 2
/// @brief Orientation distribution. Align all particles along Z axis.
# define ALIGNZ 3
/// @brief Orientation distribution. Align all particles at 45 degree angle.
# define ALIGN45 4
/// @brief Orientation distribution. Align all particles randomly within XY axis.
# define PLANEZ 5
/// @brief Orientation distribution. Align all particles randomly within XZ axis.
# define PLANEY 6
/// @brief Orientation distribution. Align all particles randomly within YZ axis.
# define PLANEX 7
/// @brief Orientation distribution. Align all particles in the direction of origin towards positive right hand corner of any cartesian plane.
# define ALIGNTR 8
/// @brief Orientation distribution. Align all particles in the direction of two oppositely charged defects.
# define ALIGNDEFECTPAIR 9

/// @brief Lower cutoff for Monte Carlo method to generate new orientations from Maier-Saupe.
# define BUSMIN 0.5
/// @brief Upper cutoff for Monte Carlo method to generate new orientations from Maier-Saupe.
# define BUSMAX 5.0

/* ****************************************** */
/* *************** NEMATIC LC *************** */
/* ****************************************** */
/// @brief Liquid crystal setting. Isotropic fluid --- not a liquid crystal.
# define ISOF 0
/// @brief Liquid crystal setting. Nematic LC using the local S value.
# define LCL 1
/// @brief Liquid crystal setting. Nematic LC using the global S value.
# define LCG 2
/// @brief Value for bacterial simulations. Make sure particles do not rotate too much.
# define BCT 3

/* ****************************************** */
/* ******** HYDRODYNAMIC INTERACTIONS ******* */
/* ****************************************** */
/// @brief Hydrodynamic interactions turned off.
# define HIOFF 1
/// @brief Hydrodynamic interactions left on.
# define HION 0

/* ****************************************** */
/* ************** NON-IDEAL EOS ************* */
/* ****************************************** */
/// @brief Non-ideal gas EOS settings. Ideal gas equation of state.
# define INCOMPOFF 0
/// @brief Non-ideal gas EOS settings. Incompressibility by swapping velocities.
# define INCOMPSWAP 1
/// @brief Non-ideal gas EOS settings. Incompressibility by virial coefficients on the local density.
# define INCOMPVIRIAL 2
/// @brief Non-ideal gas EOS settings. Incompressibility by subtracting the calculated divergence.
# define INCOMPSUB 3

/* ****************************************** */
/* ********* MULTIPHASE INTERACTIONS ********* */
/* ****************************************** */
/// @brief Multiphase interaction settings. Multiphase interactions left off.
# define MPHOFF 0
/// @brief Multiphase interaction settings. Kira Koch's Surface fitter-version of multiphase interactions.
# define MPHSURF 1
/// @brief Multiphase interaction settings. Point gradient-version of multiphase interactions.
# define MPHPOINT 2

/* ****************************************** */
/* *************** THERMOSTAT *************** */
/* ****************************************** */
//The following global variables are flags for the thermostat used
//Used by TSTECH in inputList in c-code and 'tsTech' in json input
/// @brief Thermostat settings. NOTHERM indicates that no thermostat is implemented.
# define NOTHERM 0
/// @brief Thermostat settings. VSC indicates that velocity scaling is used as the thermometer.
# define VSC 1
/// @brief Thermostat settings. BEREND indicates that the Berendsen thermostat is used.
# define BEREND 2
/// @brief Thermostat settings. HEYES indicates that Heyes thermostat for cell by cell thermostatting used.
# define HEYES 3
/// @brief Thermostat settings. MAXV is not really a thermostat. Use it with RTECH=MPCAT to have a maximum velocity vector (which is inputted as GRAV[] in input.inp).
# define MAXV 4

/* ****************************************** */
/* ********** POLYMER TRANSLOCATION ********* */
/* ****************************************** */
//For translocation, the pore width can be changed in the BC inputs
//HOWEVER, for initializing the polymer conformation, we need to know the pore width and it is hardcoded as transPoreWidth
# define transPoreWidth 2

/* ****************************************** */
/* *********** COLLISION OPERATOR *********** */
/* ****************************************** */
//The following global variables are flag values for the collision operator / rotation technique.
//Used by RTECH in inputList in c-code and 'collOp' (or 'rTech') in json input
/// @brief Collision operator option. ARBAXIS indicates that the rotation operator is a rotation about a randomly chosen axis.
# define ARBAXIS 0
/// @brief Collision operator option. ORTHAXIS applies the rotatation about one of the three cartesian axes (randomly chosen with the sign of the rotation angle also random.
# define ORTHAXIS 1
/// @brief Collision operator option. AT indicates that the Andersen Thermostat version of MPC is used.
# define MPCAT 2
/// @brief Collision operator option.  The Andersen version that conserves angular momentum.
# define RAT 3
/// @brief Collision operator option.  The Langevin version of MPCD.
# define LANG 4
// Killed the two specific no HI collision rules and the multiphase collision rule. 
// Now made these general for any collision rule. They can now be used overlayed on ANY collision operator. 
// 5, 6 and 18 are now available to be re-used as new collision operators.
// # define NOHI_ARBAXIS 5
// # define NOHI_MPCAT 6
// # define XXXX_MULTIPHASE_XXXX 18
/// @brief Collision operator option. An MPCD version of the Vicsek algorithm.
# define VICSEK 7
/// @brief Collision operator option. An MPCD version of the Chate algorithm.
# define CHATE 8
/// @brief Collision operator option. Active SRD with ARBAXIS.
# define ACT_ARBAXIS 9
/// @brief Collision operator option. Active SRD with ORTHAXIS.
# define ACT_ORTHAXIS 10
/// @brief Collision operator option. Standard Vicsek active Andersen Thermostat version of MPCD.
# define VICSEK_MPCAT 11
/// @brief Collision operator option. Standard Vicsek active Langevin version of MPCD.
# define VICSEK_LANG 12
/// @brief Collision operator option. Chate active nematic Andersen Thermostat version of MPCD.
# define CHATE_MPCAT 13
/// @brief Collision operator option. Chate active nematic Langevin version of MPCD.
# define CHATE_LANG 14
/// @brief Collision operator option. Cell-based dipole force in direction of centre of mass velocity.
# define DIPOLE_VCM 15
/// @brief Collision operator option. Cell-based dipole force in direction of local director (sum of all activities).
# define DIPOLE_DIR_SUM 16
/// @brief Collision operator option. Cell-based dipole force in direction of local director (average of all activities).
# define DIPOLE_DIR_AV 17
/// @brief Collision operator option. The Langevin version of MPC that conserves angular momentum.
# define RLANG 19
/// @brief Collision operator option. Cell-based dipole force in direction of local director (average of all activities with sigmoidal falloff).
# define DIPOLE_DIR_SIG 20
/// @brief Collision operator option. Cell-based dipole force in direction of local director (sum of all activities with sigmoidal falloff).
# define DIPOLE_DIR_SIG_SUM 21

/* ****************************************** */
/* ******************** BCS ***************** */
/* ****************************************** */
//The following global variables are different boundary conditions: COLL_TYPE
/// @brief BC type. Respects all conservation laws by impact analysis.
# define BC_IMP 0
/// @brief BC type. Rule based collisions at the suface.
# define BC_SURF_RULES 1
/// @brief BC type. Probabilistic reflections at the surface.
# define BC_THERMO_SURF 2
/// @brief BC type. Rule based collisions at half the time step.
# define BC_HALF_RULES 3
/// @brief BC type. Probabilistic reflections at half the time step.
# define BC_THERMO_HALF 4
// /// @brief BC type. Moving wall.
// #define BC_MOVING_WALL 5
/// @brief BC surface type option for SURF_MODE. The TRAD_SURF is the smooth surface of the form
/// - \f$ \mathcal{S} = \left( \frac{x-Q_x}{A_x} \right)^{p_x} + \left( \frac{y-Q_y}{A_y} \right)^{p_y} + \left( \frac{z-Q_z}{A_z} \right)^{p_z} - R^{p_R} = 0 
# define TRAD_SURF 0
/// @brief BC surface type option for SURF_MODE. The VERT_SURF is Alex Houston's vertices-based approach
# define VERT_SURF 1

/* ****************************************** */
/* ********** Apply BC to particle ********** */
/* ****************************************** */
/// @brief BC interactions left on.
# define BCON 1

/* ****************************************** */
/* ************** MD Coupling *************** */
/* ****************************************** */
/// @brief MD coupling option. No MD coupling.
# define noMD 0
/// @brief MD coupling option. MD particles included in MPCD collision operator.
# define MDinMPC 1
/// @brief MD coupling option. MPCD particles are treated in pair interactions with MD particles, but not between each other.
# define MPCinMD 2

/* ****************************************** */
/* ************** MD Warmup ***************** */
/* ****************************************** */
/// @brief MD warmup option. No MD integration and no coupling to MPCD during MPCD warmup.
# define FROZEN_WARMUP 0
/// @brief MD warmup option. MD coupled and integrated during MPCD warmup.
# define FREE_WARMUP 1
/// @brief MD warmup option. MD coupled and integrated during MPCD warmup but its center is kept where it is initiated. 
# define PINNED_WARMUP 2

/* ****************************************** */
/* ************** Monte Carlo *************** */
/* ****************************************** */
// annealNum=(int)(MCINT+MCSLOPE*effM*S/KBT);
/// @brief Monte Carlo setting. Slope of Monte Carlo.
# define MCSLOPE 0.1
/// @brief Monte Carlo setting. Interval between Monte Carlo steps.
# define MCINT 5
/// @brief Number of pseudo-particles per MPCD used to calculate the accessible volume by Monte Carlo integration. 
# define NUMMC 1000

/* ****************************************** */
/* **************** Swimmers **************** */
/* ****************************************** */
/// @brief Swimmer type. Swimmer fixed in place.
# define DUMBBELL_FIXED 0
/// @brief Swimmer type. Moving swimmers with no excluded volume interactions.
# define DUMBBELL 1
/// @brief Swimmer type. Excluded volume interactions between pairs of swimmers.
# define DUMBBELL_EXVOL 2
/// @brief Swimmer type. Only one force on the head, in the direction of the head, producing a monopole force.
# define DUMBBELL_MONOF 3
/// @brief Swimmer type. Swimmers constrained to be near the y=0 (2D) or z=0 (3D) domain boundary.
# define DUMBBELL_NEARWALL 4
/// @brief Swimmer type. Model of a <i>C. elegan</i> worm.
# define UNDULATOR 5
//If the swimmer is in the run or tumble phase
/// @brief Swimmer run-tumble phase. Running.
# define RUNNING 0
/// @brief Swimmer run-tumble phase. Tumbling.
# define TUMBLING 1
/// @brief Swimmer run-tumble phase. Shrinking.
# define SHRINKING 2
/// @brief Swimmer run-tumble phase. Extending.
# define EXTENDING 3
//Dumbbell spring
/// @brief Swimmer spring type. FENE spring.
# define FENESPRING 0
/// @brief Swimmer spring type. Hookes spring.
# define HOOKESPRING 1

/* ****************************************** */
/* ******************* DEBUG **************** */
/* ****************************************** */
//The following global variables are flag values for debugging
//Used by DBUG in c-code and 'debugOut' in json input
//Debug modes:
//Often used:
/// @brief Debug/verbosity level. Outputs nothing but sim start and sim end.
# define DBGRUN 0
/// @brief Debug/verbosity level. Only outputs warnings.
# define DBGWARN 1
/// @brief Debug/verbosity level. Only outputs initialisation information.
# define DBGINIT 2
/// @brief Debug/verbosity level. Only outputs iteration steps.
# define DBGSTEPS 3
/// @brief Debug/verbosity level. Gives a title to some methods and outputs them when called.
# define DBGTITLE 4
/// @brief Debug/verbosity level.
# define DBGWAIT 5
/// @brief Debug/verbosity level.
# define DBGTHERM 6
/// @brief Debug/verbosity level.
# define DBGHIST 7
/// @brief Debug/verbosity level.
# define DBGBCCNT 8
/// @brief Debug/verbosity level.

//Rarely used:
/// @brief Debug/verbosity level.
# define DBGMPCBC 9
/// @brief Debug/verbosity level.
# define DBGBCMPC 10
/// @brief Debug/verbosity level.
# define DBGBCBC 11
/// @brief Debug/verbosity level.
# define DBGBCMAX 12
/// @brief Debug/verbosity level.
# define DBGLCCOL 13
/// @brief Debug/verbosity level.
# define DBGBCORI 14
/// @brief Debug/verbosity level.
# define DBGJEFF 15
/// @brief Debug/verbosity level.
# define DBGMAG 16
/// @brief Debug/verbosity level.
# define DBGBINARY 17
/// @brief Debug/verbosity level.
# define DBGSWIMMER 18
/// @brief Debug/verbosity level.
# define DBGRUNTUMBLE 19
/// @brief Debug/verbosity level.
# define DBGESCAPE 20
/// @brief Debug/verbosity level.
# define DBGSWIMMERDEETS 21
/// @brief Debug/verbosity level.
# define DBGSWIMMERTORQUE 22
/// @brief Debug/verbosity level.
# define DBGINCOMP 23
#endif
