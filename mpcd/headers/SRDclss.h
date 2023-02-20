///
/// @file
///
/// @brief Contains all the structs (pseudo-classes) used in the code.
///
/// The code makes heavy use of C typedef'd structs. This file contains all of the structs
/// used throughout the code.
///
/// Members of these structs usually have a 1:1 correspondence with JSON input file parameters, and are displayed when
/// appropriate.
///

#ifndef SRDCLSS_H
#define SRDCLSS_H

# include "definitions.h"
# include "../../md/mdtypes.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ Program Classes ************* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
///
/// @brief A struct representing a single MPCD particle.
///
/// The central structure used throughout MPCD. 
/// Particles are stored in a single array.
/// However, they are commonly accessed through linked lists belonging to each cell, and then
/// pointed to by all other structures.
///
typedef struct particleMPC {
	double Q[3];			        ///< The particle position vector.
	double V[3];			        ///< The particle velocity vector.
	double U[3];			        ///< The particle alignment/orientation vector.
	double T[3];			        ///< Torque due to shear and magnetic field.
	int S_flag;				        ///< Integer flag for streaming status. 0 if must stream. 1 if has streamed.
	int SPID;					    ///< The species ID in the spec array, corresponding to this particle's species.
	double q;					    ///< Charge of  MPCD particle (for hybrid MD).
	struct particleMPC *next;	    ///< Pointer to next particle in the cell list.
	struct particleMPC *previous;	///< Pointer to previous particle in the cell list.
} particleMPC;

///
/// @brief A struct representing the hyper-parameters for a single species of MPCD fluid particles.
///
/// Contains common parameters for all particles of a given MPCD fluid species (or type or particle). 
/// These are stored in an array of these structures.
/// The these are set through json input files and read in via readJson().
/// @see readJson()
///
typedef struct spec {
	int POP;	            ///< Total population of the MPCD particles --- json `'pop'` (`'dens'` can override it for simple geometries).
	int QDIST;				///< How the position of this species' particles is initialised --- json `'qDist'`.
	int VDIST;				///< How the velocity of this species' particles is initialised --- json `'vDist'`.
	int ODIST;				///< How the orientation of this species' particles is initialised --- json `'oDist'`.
	double MASS;            ///< Mass of this species' particles --- json `'mass'`.
	double RFC;				///< Rotational friction coefficient of species nematogens --- json `'rfc'`.
	double TUMBLE;		    ///< Tumbling parameter. abs(TUMBLE)<1 for tumbling LC, abs(TUMBLE)>1 for flow aligning --- json `'tumble'`.
				            		// Can be related to effective aspect ratio of the nematogens (p) by TUMBLE=(p^2-1)/(p^2+1) 
	double CHIHI;			///< A susceptibility to shear. Not theoretical; just practical --- json `'shearSusc'`.
	double CHIA;			///< Magnetic susceptibility anisotropy chi_parallel-chi_perpendicular --- json `'magnSusc'`.
	double LEN;				///< Effective rod length to couple torque on  MPCD into force on BC (smaller=>stronger; bigger=>weaker) --- json `'len'`.
	double ACT;				///< The activity of the particles --- json `'act'`.
	double DAMP;			///< A damping/friction coefficient to go from wet to dry (to kill hydrodynamics) [0,1] --- json `'damp'`.
	double M[MAXSPECI];	    ///< Interaction matrix for multiphase fluids. Each species has a different interaction with all others --- json `'interMatr'`.
	double SIGWIDTH;		///< The width of the sigmoid for active dipole sigmoid (`DIPOLE_DIR_SIG` in definitions.h). 
	double SIGPOS;			///< The position of the sigmoid for active dipole sigmoid (`DIPOLE_DIR_SIG` in definitions.h).
	double MINACTRATIO;		///< Minimum proportion of particles in the cell for activity to be applied. 
} spec;

///
/// @brief A struct representing a single simulation boundary/wall.
///
/// Contains the information required for BC calculations. 
/// This defines <b>both</b> the boundary surface and the rules for particle interactions with the boundary conditions. 
/// All the walls are stored as an array of these structures. 
/// A description of how boundary surfaces are calculated from these prameters is provided. 
/// The these are set through json input files and read in via readJson(). 
/// @see readJson()
///
/// @note
/// Setting fixed boundary surfaces or initializing mobile surfaces uses the variables:
/// - `Q[3]` for position, \f$\vec{Q}\f$. 
/// - `AINV[3]` the "semi-axes" of the surface (as in an ellipsoid), which are the coefficients in front of each cartesion term, \f$\vec{A}\f$. 
/// - `R` for the shift term (as in the radius of a sphere), \f$ R \f$. 
/// - `P[4]` the powers  of each cartesion term, \f$\vec{p}\f$. 
///

/// @note
/// Boundary surface take on a general form of
/// - \f$ \mathcal{S} = \left( \frac{x-Q_x}{A_x} \right)^{p_x} + \left( \frac{y-Q_y}{A_y} \right)^{p_y} + \left( \frac{z-Q_z}{A_z} \right)^{p_z} - R^{p_R} = 0 \f$ where the 
/// 	+ position of the wall is set by \f$ \vec{Q}=[Q_x,Q_y,Q_z] \f$, which has the variable name `Q[3]`
/// 	+ semiaxes \f$ \vec{A}=[A_x,A_y,A_z] \f$, which has the variable name `AINV[3]`
/// 	+ radius (or scalar shift) \f$ R \f$, which has the variable name `R`
/// 	+ powers on each cartesian term \f$ \vec{p} \f$, which has the variable name `P`
/// - The surface equation becomes `( A[0]*(x-Q[0]) )^P[0] + (A[1]*(y-Q[1]))^P[1] + (A[2]*(z-Q[2]))^P[2] - R^P[3] = 0`.
///
/// Some common examples:
///
/// - PLANE:    `A=normal`, `p=[1,1,1,1]`, `Q=point on plane`, `R=0`
/// - CIRCLE:   `A=[1,1,1]`, `p=[2,2,2,2]`, `Q=centre`, `R=radius`
/// - ELLIPSE:  `A=focii`, `p=[2,2,2,2]`, `Q=centre`, `R=1`
/// - SQUIRCLE: `A=[1,1,1,1]`, `p>=4` (even), `Q=centre`, `R=radius`
///
///	There are additional complications:
///	1. An absolute operator can be used by turning on the `ABS` flag, which puts absolute operators around the terms:
///		- \f$ \mathcal{S} = \left| \frac{x-Q_x}{A_x} \right|^{p_x} + \left| \frac{y-Q_y}{A_y} \right|^{p_y} + \left| \frac{z-Q_z}{A_z} \right|^{p_z} - R^{p_R} = 0 \f$
///	2. This only allows a certain subset of rotational symmetries (even and less than or equal to 4-fold). For
///     different symmetries, we have ROTSYMM[2]. Let:
///     - `r=sqrt( (x-Q[0])^2 + (y-Q[1])^2 + (z-Q[2])^2 )`
///		- `phi=atan2( (y-Q[1])/(x-Q[0]) )`
///		- `theta=arccos( (z-Q[2])/r )`
///		- `abs( A[0]*cos(ROTSYMM[0]*phi/4)*sin(ROTSYMM[1]*theta/4) )^P[0]`
///				 `+ abs( A[1]*cos(ROTSYMM[0]*phi/4)*sin(ROTSYMM[1]*theta/4) )^P[1]`
///					`+ abs( A[2]*cos(ROTSYMM[0]*phi/4)*sin(ROTSYMM[1]*theta/4) )^P[2] - (R/r)^P[4] = 0`
///
///	This formulation allows us to make things like triangles, hexagons, stars.
/// Please see sampleInputs/ for more examples. 
///
typedef struct bc {
	// Variables that set the geometry of the surface
	double P[4];			///< Boundary `P` parameter --- json `'P'`.
	double A[3];	        ///< Boundary 'A' parameter.
    double AINV[3];         ///< Inverse of boundary 'A' parameter --- json `'aInv'` (aInv is prefered for ellipse and similar surfaces).
	double B[3];			///< Generalized amplitudes and frequencies for wavy-walls --- json `'wavy'`.
	double R;				///< Boundary 'R' parameter --- json `'R'`.
	int PLANAR;				///< Flags if BC is just a simple (X,Y or Z) plane.
	int REORIENT;			///< Flags whether or not a rotation needs to be done everytime (very expensive).
	int ABS;				///< Flags if each term should be magnitude only --- json `'abs'`.
	double ROTSYMM[2];		///< Sets the rotational symmetry of the shapes --- json `'rotSym'`.
	                        // (see Gielis, American Journal of Botany 90(3): 333–338. 2003)
	int INV;				///< Flags whether the BC is inversed or not (0-no; 1-yes) --- json `'inv'`.

	// Variables that control collisions with the wall
	//N stands for normal to surface, T for tangential
	int COLL_TYPE;		    ///< See the list of different types of collisions --- json `'colType'`.
	int PHANTOM;		    ///< The flag for using phantom particles at walls. Use phantom particles if 1; don't if 0 --- json `'phantom'`.
	double E;			    ///< Coefficient of restitution --- json `'E'`.
	double DN;              ///< Normal shift of a particle's position if it passes the bc --- json `'DN'`.
    double DT;		        ///< Tangential shift of a particle's position if it passes the bc --- json `'DT'`.
	double DVN;             ///< Normal shift of a particle's velocity if it passes the bc --- json `'DVN'`.
    double DVT;		        ///< Tangential shift of a particle's velocity if it passes the bc --- json `'DVT'`.
	double DVxyz[3];	    ///< The amount a particle's velocity is shifted if it passes the bc in cartesian coordinates --- json `'DVxyz'`.
	double MVN;             ///< Normal multiplier of a particle's velocity if it passes the bc --- json `'MVN'`.
    double MVT;		        ///< Tangential multiplier of a particle's velocity if it passes the bc --- json `'MVT'`.
	double MUN;             ///< Normal multiplier of a particle's orientation if it passes the bc --- json `'MUN'`.
    double MUT;		        ///< Tangential multiplier of a particle's orientation if it passes the bc --- json `'MUT'`.
	double MUxyz[3];	    ///< The amount a particle's orientation is multiplied by if it pass the bc --- json `'MUxyz'`.
	double DUxyz[3];	    ///< The amount added to a particle's orienation if it passes the bc in cartesian coordinates --- json `'DUxyz'`.
	double KBT;			    ///< The temperature of the wall (only used if COLL_TYPE is set to thermal collisions) --- json `'KBT'`.

	// The surfaces coordinates
	double Q[3];		    ///< The BC's position --- json `'Q'`.
	double V[3];		    ///< The BC's velocity --- json `'V'`.
	double O[3];		    ///< The BC's orientation (angle about x,y,z) --- json `'O'`.
	double L[3];		    ///< The BC's angular velocity --- json `'L'`.
	double G[3];		    ///< The BC's external acceleration --- json `'G'`.
	double W;			    ///< The particle's W for passing this boundary.
	double Q_old[3];	    ///< The BC's last position.
	double O_old[3];	    ///< The BC's last orientation.
	double dV[3];		    ///< The BC's change in velocity due to collisions.
	double dL[3];		    ///< The BC's change in angular velocity due to collisions.

	// Qualities of the object
	int DSPLC;				///< Flags whether BC can be pushed around by particles (0-no; 1-yes) --- json `'dsplc'`.
	double MASS;			///< The BC's mass (only relevant if it moves) --- json `'mass'`.
	double VOL;				///< Body's volume.
	double I[3][3];		    ///< The body's moment of inertia.
/*
   Examples

   A periodic boundary condition: Top wall
   COLL_TYPE = {0,1}		Either impulse or rule based
   PHANTOM = 1			Don't want film of less viscous fluid surrounding CV
   E = -1.			*** Since the particle it not transfering momentum
   Q = (0,0,0)			Plane's origin
   V = (0,0,0)			At rest
   L = (0,0,0)			At rest
   A = (0,-1,0)			Normal must point INTO the control volume
   P = 1			Planar
   R = 150			Plane's distance from Q along A
   DN = 150, DT = 0		Shift position 150 along normal, don't shift tangential to plane
   DVT = 0, DVN = 0		Don't shift velocity
   DVxyz =[0,0,0]		Don't shift velocity
   MVT = 1, MVN = 1		Don't change direction
   INV = 0			Don't inverse the normal.
   M = ANYTHING
   KBT = ANYTHING
		Putting these into the equation we have y = 150.

   2D Sphere: radius 2
   COLL_TYPE = ANYTHING
   PHANTOM = ANYTHING
   E = 1.			Elastic
   Q = (75,75,0)		Centre of a 150X150X0 control volume
   V = (0,0,0)			Starts at rest
   L = (0,0,0)			Starts at rest
   A = (1,1,0)			A[2] = 0 cuz 3D
   P = 2			Sphere/ellipse
   R = 4			Square of radius og 2
   DN = 0, DT = 0		Don't shift position
   DVT = 0, DVN = 0		Don't shift velocity
   DVxyz =[0,0,0]		Don't shift velocity
   MVT = {1,-1}, MVN = {-1,-1}	For reflective use {MTV=0,MVN=-1} so only normal flips and tang unchanged. For bounce back use {MVT=-1,MVN=-1}
   INV = 0			Don't inverse the normal.
   M = ANYTHING
   KBT = ANYTHING
		Putting these into the equation we have (x-75)^2 + (y-75)^2 = 4

   Hard wall: left wall
   COLL_TYPE = ANYTHING
   PHANTOM = 1			Don't want slip
   E = 1.			Elastic
   Q = (0,0,0)			Plane's origin
   V = (0,0,0)			At rest
   L = (0,0,0)			At rest
   A = (1,0,0)			Normal must point INTO the control volume
   P = 1			Planar
   R = 0			Plane's distance from Q along A
   DN = 0, DT = 0		Shift position 150 along normal, don't shift tangential to plane
   DVT = 0, DVN = 0		Don't shift velocity
   DVxyz =[0,0,0]		Don't shift velocity
   MVT = {1,-1}, MVN = {-1,-1}	Reflect or bounce back
   INV = 0			Don't inverse the normal.
   M = ANYTHING
   KBT = ANYTHING
		Putting these into the equation we have x = 0.

   Capillary: radius 10
   COLL_TYPE = ANYTHING
   PHANTOM = 1			Don't want slip
   E = 1.			Elastic
   Q = (0,Y_CENTRE,Z_CENTRE)	Plane's origin
   V = (0,0,0)			At rest
   L = (0,0,0)			At rest
   A = (0,1,1)		Normal must point INTO the control volume
   P = 2			Planar
   R = RADIUS^2=100			radius squared
   DN = 0, DT = 0		Shift position 150 along normal, don't shift tangential to plane
   DVT = 0, DVN = 0		Don't shift velocity
   DVxyz =[0,0,0]		Don't shift velocity
   MVT = {1,-1}, MVN = {-1,-1}	Reflect or bounce back
   INV = 1			DO inverse the normal so that the fluid particles are INSIDE.
   M = ANYTHING
   KBT = ANYTHING
		If INV = 0 then we have a pillar/rod rather than a capillary

*/
} bc;

///
/// @brief The struct representing a single MPCD cell
///
/// This struct acts as a container for all information about a single MPCD cell. These are primarily stored as a 3D
/// array representing the system domain.
///
typedef struct cell {
	int POP;					///< Total population of the cell (MPCD particles, swimmers, and monomers).
	int POPSRD;					///< Population of MPCD particles.
	int POPSW;					///< Population of swimmers.
	int POPMD;					///< Population of MD particles.
	double MASS;				///< Total mass of the cell.
	int SP[MAXSPECI];			///< Subpopulations of each MPCD species type in the cell.
	double I[3][3];		        ///< The cell's moment of inertia about (0,0,0).
	double E[3][3];		        ///< Velocity gradient tensor. First index `[i]` is on velocity, second `[j]` on derivative: `E[i][j]= dv[i]/dx[j]`.
	double S;					///< The cell's scalar order parameter.
	double CM[3];			    ///< Centre of mass position of the cell.
	double VCM[3];		        ///< Centre of mass velocity of the cell.
	double DIR[3];		        ///< Director of the cell (average orientation).
	double FLOW[3];		        ///< Centre of mass velocity of the cell averaged over FLOWOUT time steps.
	double Ps[3][3];		    ///< Streaming part of the local instantaneous stress tensor.
	double Pc[3][3];		    ///< Collisional part of the local instantaneous stress tensor.

	struct particleMPC *pp;		///< Pointer to first SRD particle in list.
	struct particleMD *MDpp;	///< Pointer to first MD particle in list.
	struct smono *sp;			///< Pointer to first swimmer monomer particle in list.
} cell;

///
/// @brief Helper container struct containing all output files.
///
/// All output files are stored within this struct, used as a container for simple passing to functions.
///
typedef struct outputFilesList {
	FILE *fcoarse,*fflow,*fenergy,*fenergyfield,*fenneighbours;
	FILE *fsynopsis,*favvel,*forder,*forderQ,*forderQK,*favs,*fdensSTD,*fchckpnt,*fenstrophy,*fmultiphase,*fpressure;
	FILE *fcorrVV,*fcorrNN,*fcorrWW,*fcorrDD,*fcorrSS,*fcorrPP,*fbinder;
	FILE *fhistVel,*fhistSpeed,*fhistVort,*fhistEnstr,*fhistDir,*fhistS,*fhistDens;
	FILE *fenergyspect,*fenstrophyspect;
	FILE *ftopo,*fdefects,*fdisclination;
	FILE *fdetail[MAXSPECI];
	FILE *fsolids[MAXBC];
	FILE *fswimmers,*fswimmersOri,*fruntumble;
} outputFilesList;

///
/// @brief Container for possible outputs. Each int represents the time period between dumps of that output.
///
/// This struct is used to store the time period between dumps of each output. Each dump is stored as an integer: A
/// value of 0 means that no dump takes place, a finite positive value indicates the period between dumps for that
/// output.
/// The these are set through json input files and read in via readJson().
/// @see readJson()
///
typedef struct outputFlagsList {
	int TRAJOUT;				///< Flag for if the detailed trajectories of every particle are outputted --- json `'trajOut'`.
	int COAROUT;				///< Flag for if coarse grain is outputted --- json `'coarseOut'`.
	int AVVELOUT;				///< Flag for if total average velocity is outputted --- json `'avVelOut'`.
	int FLOWOUT;				///< Flag for if the flow field is outputted --- json `'flowOut'`.
	int HISTVELOUT;             ///< Flag for if the velocity distribution is outputted --- json `'histVelOut'`.
    int HISTSPEEDOUT;           ///< Flag for if the speed distribution is outputted --- json `'histSpeedOut'`.
    int HISTVORTOUT;            ///< Flag for if the vorticity distribution is outputted --- json `'histVortOut'`.
    int HISTENSTROUT;           ///< Flag for if the enstrophy distribution is outputted --- json `'histEnsOut'`.
    int HISTDIROUT;             ///< Flag for if the director distribution is outputted --- json `'histDirOut'`.
    int HISTSOUT;               ///< Flag for if the scalar order distribution is outputted --- json `'histSOut'`.
    int HISTNOUT;	            ///< Flag for if the number per cell distribution is outputted --- json `'histNOut'`.
	int ENERGYSPECTOUT;         ///< Flag for if energy spectra are outputted --- json `'energySpecOut'`.
    int ENSTROPHYSPECTOUT;	    ///< Flag for if enstrophy spectra are outputted --- json `'enstrophySpecOut'`.
	int TOPOOUT;                ///< Flag for if topological charge field is outputted --- json `'topoFieldOut'`.
    int DEFECTOUT;              ///< Flag for if defect position list is outputted --- json `'defectsOut'`.
    int DISCLINOUT;	            ///< Flag for if disclination tensor field are outputted --- json `'disclinOut'`.
	int ENOUT;					///< Flag for if system energy is outputted --- json `'energyOut'`.
	int ENFIELDOUT;             ///< Flag for if orientational energy field is outputted --- json `'oriEnOut'`.
    int ENNEIGHBOURS;	        ///< Flag for if orientational energy as a function of neighbours is outputted --- json `'neighbourEnOut'`.
	int SPOUT;					///< Flag for if the colour/phi/species-type field is outputted --- json `'colourOut'`.
	int PRESOUT;				///< Flag for if the pressure field is outputted --- json `'pressureOut'`.
	int SYNOUT;					///< Flag for if system synopsis is outputted --- json `'synopsisOut'`.
	int ORDEROUT;				///< Flag for if the order parameter are outputted --- json `'dirSOut'`.
	int QTENSOUT;               ///< Flag for if the order parameter tensor is outputted --- json `'qTensOut'`.
    int QKOUT;			        ///< Flag for if the reciprocal space Q is outputted --- json `'qkTensOut'`.
	int AVSOUT;					///< Flag for if total average scalar order parameter is outputted --- json `'avSOut'`.
	int DENSOUT;				///< Flag for density standard deviation is outputted --- json `'densSDOut'`.
	int ENSTROPHYOUT;			///< Flag for if total average enstrophy is outputted --- json `'enstrophyOut'`.
	int CHCKPNT;				///< Flag for checkpointing --- json `'checkpointOut'`.
	float CHCKPNTTIMER;			///< Flag for checkpointing based on a timer --- json `'checkpointTimerOut'`.
	int CHCKPNTrcvr;			///< Flag for simulation from recovery of checkpoint.
	int BINDER;                 ///< Flag for Binder cumulant --- json `'binderOut'`.
    int BINDERBIN;		        ///< Flag for the bin size of the Binder cumulant --- json `'binderBin'`.
	int printSP;				///< How many of the species are printed --- json `'trajSpecOut'`.
	int SOLOUT;					///< How often the solid bc trajectories are outputted --- json `'solidTrajOut'`.
	int CVVOUT;                 ///< Flag for if velocity-velocity spatial correlations are outputted --- json `'velCorrOut'`.
    int CNNOUT;                 ///< Flag for if director-director spatial correlations are outputted --- json `'dirCorrOut'`.
    int CWWOUT;	                ///< Flag for if vorticity-vorticity spatial correlations are outputted --- json `'vortCorrOut'`.
	int CDDOUT;                 ///< Flag for if density-density spatial correlations are outputted --- json `'densCorrOut'`.
    int CSSOUT;                 ///< Flag for if order-order spatial correlations are outputted --- json `'orderCorrOut'`.
    int CPPOUT;	                ///< Flag for if phi-phi (phase) spatial correlations are outputted --- json `'phaseCorrOut'`.
	int SWOUT;                  ///< Flag for if swimmer positions are outputted --- json `'swimQOut'`.
    int SWORIOUT;               ///< Flag for if swimmer orientations are outputted --- json `'swimOOut'`.
    int RTOUT;	                ///< Flag for if swimmer run/tumbles are outputted --- json `'swimRTOut'`.
} outputFlagsList;

///
/// @brief Helper container struct that contains kinetic theory parameters.
///
/// This container structure is used to store all fluid kinetic theory parameters that may be of use, to make it easy
/// to pass them around.
///
typedef struct kinTheory {
	double MFP;					///< Mean Free Path.
	double VISC;				///< Kinematic viscosity.
	double THERMD;				///< Thermal Diffusion.
	double SDIFF;				///< Self diffusion.
	double SPEEDOFSOUND;		///< Speed of sound.
	double sumM;				///< Sum of all masses.
} kinTheory;

///
/// @brief Helper container struct that contains input parameters.
///
/// This container structure was used to store the parameters associated with the legacy input.inp file. Now it stores
/// all of the "central" parameters of the system simulation. It is used to simplify passing these parameters around
/// between methods.
/// The these are set through json input files and read in via readJson().
/// @see readJson()
///
typedef struct inputList {
	double KBT;					///< Temperature: A third of thermal energy. Sets energy scale --- json `'kbt'`.
	double dt;					///< MPCD time step value --- json `'dt'`.
    double tolD;				///< Tolerance of defect tracker --- json `'tolD'`.
	int stepsMD;				///< Number of MD steps per SRD step (NOTE make variable) --- json `'stepsMD'`.
	int warmupSteps;            ///< Number of iterations in the warmup phase --- json `'warmUp'`.
    int simSteps;	            ///< Number of iterations in the simulation phases --- json `'simSteps'`.
	unsigned long seed;			///< Seed for random number generator --- json `'seed'`.
	double RA;                  ///< Rotation angle --- json `'rotAng'`.
    double C;                   ///< Rotation angle cos(RA).
    double S;				    ///< Rotation angle sin(RA).
	double GRAV[_3D];			///< Constant acceleration from external force --- json `'grav'`.
	double MAG[_3D];			///< Constant external magnetic field to torque nematogens --- json `'mag'`.
	int GRAV_FLAG;              ///< Flag for if no acceleration.
    int MAG_FLAG;		        ///< Flag for if no torque.
	double FRICCO;				///< Friction coefficient for Langevin thermostat --- json `'fricCoef'`.
	int TSTECH;					///< Temperature scaling technique --- json `'tsTech'`.
	double TAU;					///< The temperature relaxation time scale --- json `'tau'`.
	double MFPOT;				///< The mean-field potential from self-consistent mean-field liquid crystals --- json `'mfpot'`.
	int RTECH;					///< Rotation technique --- json `'collOp'`.
	int LC;						///< If LC=LCG=2 then liquid crystal using global S, if LC=LCL=1 then use local S, else isotropic (ISOF=0) --- json `'lc'`.
	int RFRAME;					///< Flags initial galilean trans to rest frame (0 No shift, 1 shift) --- json `'rFrame'`.
	int zeroNetMom;				///< If GAL=1 AND GRAV[D3] = [0,0,0] then zero every zeroNetMom steps. --- json `'zeroNetMom'`.
	int GALINV;					///< If rshift=1 do the random shifting of mpcd particles. If zero then do not zeroNetMom steps. --- json `'galInv'`.
	int noHI;					///< If noHI=1 remove hydrodynamic interactions by randomly scrambling of velocities.
	int inCOMP;					///< If inCOMP=1 remove div(v) by ***ALGORITHM NOT CREATED YET. BUILD ON https://doi.org/10.1063/5.0037934?***.
	int MULTIPHASE;				///< MULTIPHASE mode. If MULTIPHASE==0 then no interactions between particles of different species occurs.
} inputList;

///
/// @brief Struct that contains all of the parameters associated with a species of swimmers.
///
/// This container structure is used to store all of the hyper-parameters associated with a species of swimmers. It is
/// stored as an array.
/// The these are set through json input files and read in via readJson().
/// @see readJson()
///
typedef struct specSwimmer {
	int TYPE;					///< Type of swimmer 0=fix-dipole; 1=dumbbell; 2=dumbell with excluded vol; 3=dumbell with no counter-force  --- json `'typeSwim'`.
	int QDIST;					///< QDIST is the flag which tells how the swimmers are intially placed --- json `'qDistSwim'`.
	int ODIST;					///< ODIST is the flag which tells how the swimmers are intially placed --- json `'oDistSwim'`.
	int headM;                  ///< Mass of head monomer in the swimmer (tail mirrors the head) --- json `'headMSwim'`.
    int middM;			        ///< Mass of middle monomer in swimmer --- json `'midMSwim'`.
	double iheadM;              ///< Inverse mass of the head/tail monomers.
    double imiddM;		        ///< Inverse mass of the middle monomer.
	int HSPid;                  ///< Multiphase fluid particle type of the head monomer (tail is invisible) --- json `'hspIdSwim'`.
    int MSPid;			        ///< Multiphase fluid particle type of the middle monomer (tail is invisible) --- json `'mspIdSwim'`.
	double FS;					///< Magnitude of propulsion force to set swimming speed --- json `'fsSwim'`.
	double TS;					///< Magnitude of torque to set swimmer rotlet dipole (sign sets CW or CCW) --- json `'tsSwim'`.
	double DS;					///< Dipole strength (for a pusher DS>0; for a puller DS<0) (i.e. length of dipole in multiples of head/middle separation) --- json `'dsSwim'`.
	double sizeShrink;			///< How much Lennard Jones sigma (monomer size) and ro (equilibrium separation) are shrunk when tumbling --- json `'sizeShrinkSwim'`.
	double springShrink;		///< How much the spring constant is shrunk when tumbling --- json `'springShrinkSwim'`.
	double fixDist;				///< The fixed distance from the wall for DUMBBELL_NEARWALL mode --- json `'fixDistSwim'`.
	double k;					///< Spring constant --- json `'kSwim'`.
	double ro;                  ///< Spring separation --- json `'roSwim'`.
    double iro;				    ///< Inverse spring separation.
	double sig;                 ///< Lennard Jones sigma (monomer diameter) --- json `'sigSwim'`.
    double isig;			    ///< Inverse Lennard Jones sigma.
	double eps;					///< Lennerd Jones interaction energy --- json `'epsSwim'`.
	double runTime;             ///< Average run time in units of  MPCD time steps (dt) --- json `'runTSwim'`.
    double tumbleTime;	        ///< Average tumble time in units of  MPCD time steps (dt) --- json `'tumTSwim'`.
	int shrinkTime;				///< Set time to shrink/extend in units of  MPCD time steps (dt) --- json `'shrTSwim'`.
	double MAGMOM;				///< Magnetic moment strength/magnitude --- json `'magMomSwim'`.
} specSwimmer;

///
/// @brief Struct that represents a single monomer.
///
/// This container structure is used to store all of the parameters associated with a single monomer.
///
typedef struct smono {
	double Q[_3D];				///< Position of the monomer.
	double V[_3D];				///< Velocity of the monomer.
	double A[_3D];				///< Acceleration of the monomer.
	int HorM;					///< Signfies whether monomer is a head, or a middle (0=head; 1=middle). Tail takes the mass of the head and doesn't participate.
	struct smono *next;			///< Pointer to next particle in cell list.
	struct smono *previous;		///< Pointer to previous particle in cell list.
} smono;

///
/// @brief Struct that represents a single MPCD swimmer.
///
/// This container structure is used to store all of the parameters associated with a single MPCD swimmer.
///
typedef struct swimmer {
	smono H;                    ///< The head of the dumbbell.
    smono M;					///< The middle of the dumbbell.
	int RT;						///< Flag for if in Run phase (0) or Tumble phase (1) or Shrink phase (2) or Extension phase (3).
	double n0[_3D];				///< Orientation vector of the swimmer at start of phase.
	int timeCNT;				///< Time counter for run or tumble.
	int timeRND;				///< The randomly generated time that a given run or tumble will go for.
	double ro;                  ///< Finitely extensible nonlinear elastic (FENE) seperation --- Each swimmer needs it's own since it changes for run/tumbler.
    double iro;				    ///< Inverse FENE separation.
	double sig;                 ///< Lennard Jones size sigma --- Each swimmer needs it's own since it changes for run/tumbler.
    double isig;			    ///< Inverse Lennard Jones sigma.
	double k;					///< FENE separation --- Each swimmer needs it's own since it changes for run/tumbler.
} swimmer;

#endif
