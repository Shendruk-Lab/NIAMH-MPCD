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
// If running sims don't define DBG. If debugging define DBG and choose a debug mode (bottom of page) ---- json 'debugOut'
# define DBG
// If you want the files get printed immediately include this option for force flushing the buffer
# define FFLSH
// If you want to use the LEGACY mersenne twister RNG then uncomment below
//# define RNG_MERSENNE

/* ****************************************** */
/* ************ PROGRAM CONSTANTS *********** */
/* ****************************************** */
//The number of bins used for distributions --- best if an odd integer
# define BINS 101
//The tolerance of the precision
/* USED IN SUBROUTINES: invert2x2(), checkplacement(),checkplace(),chooseT(),secant_time() */
# define TOL 1e-8
//When doing the ACTIVE rotation in activeSRD() if the angle is too small nans result
# define ROTTOL 1E-5
//Small times have a magnitude less than SCALE
# define SCALE 1E-2
//Number of times we allow a particle to bounce between BCs in a single time step
# define NBOUNCE 50
//Maximum number of allowed MPCD species
# define MAXSPECI 10
//Maximum number of allowed MPCD species
# define MAXBC 50
//OMP chunk size for dynamic scheduling. Smaller number will result in more overhead, but less chance of thread starvation
# define OMPCHUNK 5

/* ****************************************** */
/* *********** UNIVERSAL CONSTANTS ********** */
/* ****************************************** */
# define pi M_PI
# define e M_E

/* ****************************************** */
/* ************* DIMENSIONALITY ************* */
/* ****************************************** */
# define _3D 3
# define _2D 2
# define _1D 1

/* ****************************************** */
/* ************* RANDOM NUMBER ************** */
/* ****************************************** */
# define NN 624
# define MM 397
# define MATRIX_A 0x9908b0dfUL   /* constant vector a */
# define UPPER_MASK 0x80000000UL /* most significant w-r bits */
# define LOWER_MASK 0x7fffffffUL /* least significant r bits */

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
//If set of data should be printed out to a file set flag to OUT, if not printed set to NOUT
# define OUT 1
# define NOUT 0

/* ****************************************** */
/* ***************** STREAM ***************** */
/* ****************************************** */
//If the particle needs to stream its STREAM, if not printed set to NSTREAM
# define STREAM 1
# define NO_STREAM 0

/* ****************************************** */
/* ******************* QDIST **************** */
/* ****************************************** */
//The following global variables are flag values for the variable QDIST in species structures 
//Used for both MPCD particles and swimmers --- qDist in json
//PPF indicates that the particles can be placed randomly and freely within the control volume
# define PRF 0
//READ indicates that the particles' positions should be read from a file
#define READ 1
//CENTRED places a swimmer at the centre of the control volume
#define CENTRED 2

/* ****************************************** */
/* ******************* VDIST **************** */
/* ****************************************** */
//The following global variables are flag values for the variable VDIST in species structures
//Used for both MPCD particles and swimmers --- vDist in json
//RANDVEL indicates that the particles' velosities are drawn from a have uniformly random velocity distribution
# define RANDVEL 0
//READ must be applied to both
//define READ 1
//HOMO indicates that the particles all have the same initial speed but different directions
# define HOMO 2
//HOMO indicates that the particles' initial velocity is given by a given value along each axis in either direction
# define HOMOAXIS 3
//GAUSS indicates that the particles' velocity follow a gaussian distribution
# define GAUSS 4

/* ****************************************** */
/* ******************* UDIST **************** */
/* ****************************************** */
//The following global variables are flag values for the variable UDIST in species structures
//Used for both MPCD particles and swimmers --- oDist in json
//RANDORIENT indicates that the particles' direction are drawn from a have uniformly random direction
# define RANDORIENT 0
//Align all particles along X axis
# define ALIGNX 1
//Align all particles along Y axis
# define ALIGNY 2
//Align all particles along Z axis
# define ALIGNZ 3
//Align all particles at 45 degree angle
# define ALIGN45 4
//Align all particles randomly within XY axis
# define PLANEZ 5
//Align all particles randomly within XZ axis
# define PLANEY 6
//Align all particles randomly within YZ axis
# define PLANEX 7
//Align all particles in the direction of origin towards positive right hand corner of any cartesian plane
# define ALIGNTR 8
//Align all particles in the direction of two oppositely charged defects
# define ALIGNDEFECTPAIR 9

//Lower cutoff for MC method to generate new orientations from Maier-Saupe
# define BUSMIN 0.5
//Upper cutoff for MC method to generate new orientations from Maier-Saupe
// # define BUSMAX 30.0
# define BUSMAX 5.0
//Isotropic fluid - not a liquid crystal
# define ISOF 0
//Nematic LC using the local S value
# define LCL 1
//Nematic LC using the global S value
# define LCG 2

/* ****************************************** */
/* *************** THERMOSTAT *************** */
/* ****************************************** */
//The folloing global variables are flags for the thermostat used 
//Used by TSTECH in inputList in c-code and 'tsTech' in json input
//NOTHERM indicates that no thermostat is implemented
# define NOTHERM 0
//VSC indicates that velocity scaling is used as the thermometer
# define VSC 1
//BEREND indicates that the Berendsen thermostat is used
# define BEREND 2
//HEYES indicates that Heyes thermostat for cell by cell thermostatting used
# define HEYES 3
//MAXTH is not really a thermostat. Use it with RTECH=MPCAT to have a maximum velocity vector (which is inputted as GRAV[] in input.inp)
# define MAXV 4

/* ****************************************** */
/* *********** COLLISION OPERATOR *********** */
/* ****************************************** */
//The following global variables are flag values for the collision operator / rotation technique.
//Used by RTECH in inputList in c-code and 'colOp' in json input
//ARBAXIS indicates that the rotation operator is a rotation about a randomly chosen axis
# define ARBAXIS 0
//ORTHAXIS applies the rotatation about one of the three cartesian axes (randomly chosen
//with the sign of the rotation angle also random
# define ORTHAXIS 1
//AT indicates that the Andersen Thermostat version of MPC is used.
# define MPCAT 2
// The Andersen version that conserves angular momentum
# define RAT 3
// The Langevin version of MPC
# define LANG 4
// The Brownian thermostat version (uses ARBAXIS version) i.e. no Hydrodynamic Interactions
# define NOHI_ARBAXIS 5
// The Brownian thermostat version (uses MPCAT version) i.e. no Hydrodynamic Interactions
# define NOHI_MPCAT 6
// An MPCD version of the Vicsek algorithm
# define VICSEK 7
// An MPCD version of the Chate algorithm
# define CHATE 8
// Active SRD with ARBAXIS
# define ACT_ARBAXIS 9
// Active SRD with ORTHAXIS
# define ACT_ORTHAXIS 10
// Standard Vicsek active Andersen Thermostat version of MPC
# define VICSEK_MPCAT 11
// Standard Vicsek active Langevin version of MPC
# define VICSEK_LANG 12
// Chate active nematic Andersen Thermostat version of MPC
# define CHATE_MPCAT 13
// Chate active nematic Langevin version of MPC
# define CHATE_LANG 14
// Cell-based dipole force in direction of centre of mass velocity
# define DIPOLE_VCM 15
// Cell-based dipole force in direction of local director (sum of all activities)
# define DIPOLE_DIR_SUM 16
// Cell-based dipole force in direction of local director (average of all activities)
# define DIPOLE_DIR_AV 17
// Romain and my new binary fluid
# define MULTIPHASE 18
// The Langevin version of MPC that conserves angular momentum
# define RLANG 19
// Cell-based dipole force in direction of local director (average of all activities with sigmoidal falloff)
# define DIPOLE_DIR_SIG 20
// Cell-based dipole force in direction of local director (sum of all activities with sigmoidal falloff)
# define DIPOLE_DIR_SIG_SUM 21

/* ****************************************** */
/* ******************** BCS ***************** */
/* ****************************************** */
//The following global variables are different boundary conditions: COLL_TYPE
//Respects all conservation laws by impact analysis
# define BC_IMP 0
//Rule based collisions at the suface
# define BC_SURF_RULES 1
//Probabilistic reflections at the surface
# define BC_THERMO_SURF 2
//Rule based collisions at t/2
# define BC_HALF_RULES 3
//Probabilistic reflections at t/2
# define BC_THERMO_HALF 4
//Moving wall
# define BC_MOVING_WALL 5

/* ****************************************** */
/* ************** MD Coupling *************** */
/* ****************************************** */
# define noMD 0
# define MDinMPC 1
# define MPCinMD 2

/* ****************************************** */
/* ************** Monte Carlo *************** */
/* ****************************************** */
// annealNum=(int)(MCINT+MCSLOPE*effM*S/KBT);
# define MCSLOPE 0.1
# define MCINT 5

/* ****************************************** */
/* **************** Swimmers **************** */
/* ****************************************** */
# define DUMBBELL_FIXED 0
# define DUMBBELL 1
# define DUMBBELL_EXVOL 2
# define DUMBBELL_MONOF 3
# define DUMBBELL_NEARWALL 4
# define UNDULATOR 5
//If the swimmer is in the run or tumble phase
# define RUNNING 0
# define TUMBLING 1
# define SHRINKING 2
# define EXTENDING 3
//Dumbbell spring
# define FENESPRING 0
# define HOOKESPRING 1

/* ****************************************** */
/* ******************* DEBUG **************** */
/* ****************************************** */
//The following global variables are flag values for debugging
//Debug modes:
//Often used:
# define DBGRUN 0
# define DBGWARN 1
# define DBGINIT 2
# define DBGSTEPS 3
# define DBGTITLE 4
# define DBGWAIT 5
# define DBGTHERM 6
# define DBGHIST 7
# define DBGBCCNT 8
//Rarely used:
# define DBGMPCBC 9
# define DBGBCMPC 10
# define DBGBCBC 11
# define DBGBCMAX 12
# define DBGLCCOL 13
# define DBGBCORI 14
# define DBGJEFF 15
# define DBGMAG 16
# define DBGBINARY 17
# define DBGSWIMMER 18
# define DBGRUNTUMBLE 19
# define DBGESCAPE 20
# define DBGSWIMMERDEETS 21
# define DBGSWIMMERTORQUE 22
#endif
