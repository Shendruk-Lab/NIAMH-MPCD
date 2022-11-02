# include "definitions.h"
# include "../../md/mdtypes.h"

#ifndef SRDCLSS_H
#define SRDCLSS_H
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ Program Classes ************* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
typedef struct particleMPC {
	double Q[3];			//Q is the position XYV refer to
	double V[3];			//V is the velocity
	double U[3];			//U is the alignment/orientation
	double T[3];			//T is the torque due to shear and magnetic field
	int S_flag;				//0 if must stream. 1 if has streamed.
	int SPID;					//SPID is this particle's species identification (links the particle to spec)
	double q;					//Charge of MPC particle (for hybrid MD)
	struct particleMPC *next;	//pointer to next particle in cell list
	struct particleMPC *previous;	//pointer to previous particle in cell list
} particleMPC;
typedef struct spec {
	int POP;					//POP is the population of the species
	int QDIST;				//QDIST is the flag which tells how the particles are intially placed
	int VDIST;				//VDIST is the flag which tells how the particles are intially placed
	int ODIST;				//ODIST is the flag which tells how the particles are intially oriented
	double MASS;					//M is the mass of particles of this type
	double RFC;				//RFC is the rotational friction coefficient the nematogens
	double TUMBLE;		/*TUMBLE is the tumbling parameter. abs(TUMBLE)<1 tumbling LC but abs(TUMBLE)>1 flow aligning.
				Can be related to effective aspect ratio of the nematogens (p) by TUMBLE=(p^2-1)/(p^2+1) */
	double CHIHI;			//A susceptibility to shear --- not theoretical; just practical
	double CHIA;			//Magnetic susceptibility anisotropy chi_parallel-chi_perpendicular
	double LEN;				//Effective rod length to couple torque on MPC into force on BC (smaller=>stronger; bigger=>weaker)
	double ACT;				//The activity of the particle
	double DAMP;			//A damping/friction coefficient to go from wet to dry (to kill hydrodynamics) [0,1]
	double M[MAXSPECI];	//Interaction matrix for multiphase fluids --- each species has a different interaction with all others
} spec;
typedef struct bc {
/*
   Boundary conditions take on a general form of
   ( (x-h)/a )^p + ((y-k)/b)^p + ((z-l)/c)^p = R^p
   which is in this code
   AINV[3] = (a,b,c)
   Q[3] = (h,l,l)
   ( A[0]*(x-Q[0]) )^P[0] + (A[1]*(y-Q[1]))^P[1] + (A[2]*(z-Q[2]))^P[2] - R^P[4] = 0
      PLANE    A=norm,p=1,R=-dist --- REMEMBER for a plane r must be neg!
      CIRCLE   A=[1,1,1] p=2, Q=centre, R=radius
      ELLIPSE  A=focii, p=2, Q=centre R=1
      SQUIRCLE A=1, p>=4, Q=centre, R=radius

	 There are additional complications:
	 (1) An absolute operator can be used around the terms
	 			abs( A[0]*(x-Q[0]) )^P[0] + abs(A[1]*(y-Q[1]))^P[1] + abs(A[2]*(z-Q[2]))^P[2] - R^P[4] = 0

	 (2) This only allows a certain subset of rotational symmetries (even and less than or equal to 4-fold).
	     For different symmetries, we have ROTSYMM[2]
			  Let r=sqrt( (x-Q[0])^2 + (y-Q[1])^2 + (z-Q[2])^2 )
				phi=atan2( (y-Q[1])/(x-Q[0]) )
				theta=arccos( (z-Q[2])/r )
		 	 	abs( A[0]*cos(ROTSYMM[0]*phi/4)*sin(ROTSYMM[1]*theta/4) )^P[0]
				 + abs( A[1]*cos(ROTSYMM[0]*phi/4)*sin(ROTSYMM[1]*theta/4) )^P[1]
					+ abs( A[2]*cos(ROTSYMM[0]*phi/4)*sin(ROTSYMM[1]*theta/4) )^P[2] - (R/r)^P[4] = 0
			This formulation allows us to make things like triangles, hexagons, stars
*/
	// Variables that set the geometry of the surface
	double P[4];
	double A[3],AINV[3];
	double R;
	int PLANAR;		//Flags if BC is just a simple (X,Y or Z) plane.
	int REORIENT;	//Flags whether or not a rotation needs to be done everytime (very expensive)
	int ABS;			//Flags if each term should be magnitude only
	double ROTSYMM[2];	//Sets the rotational symmetry of the shapes (see Gielis, American Journal of Botany 90(3): 333–338. 2003)
	int INV;			// Flags whether the BC is inversed or not (0-no; 1-yes)

	// Variables that control collisions with the wall
	//N stands for normal to surface, T for tangential
	int COLL_TYPE;		//See the list of different types of collisions
	int PHANTOM;			//The flag for using phantom particles at walls. Use phantom particles if 1; don't if 0
	double E;					//Coefficient of restitution
	double DN,DT;			//The amount a particle's position is shifted if it passses the bc
	double DVN,DVT;		//The amount a particle's velocity is shifted if it passes the bc
	double DVxyz[3];	//The amount a particle's velocity is shifted if it passes the bc in cartesian coordinates
	double MVN,MVT;		//The amount a particle's velocity is multiplied by if it pass the bc
	double MUN,MUT;		//The amount a particle's orientation is multiplied by if it pass the bc
	double MUxyz[3];	//The amount a particle's orientation is multiplied by if it pass the bc
	double DUxyz[3];	//The amount added to a particle's orienation if it passes the bc in cartesian coordinates
	double KBT;				//The temperature of the wall (only used if COLL_TYPE is set to thermal collisions)

	// The surfaces coordinates
	double Q[3];			//The BC's position
	double V[3];			//The BC's velocity
	double O[3];			//The BC's orientation (angle about x,y,z)
	double L[3];			//The BC's angular velocity
	double G[3];			//The BC's external acceleration
	double W;					//The particle's W for passing this boundary
	double Q_old[3];	//The BC's last position
	double O_old[3];	//The BC's last orientation
	double dV[3];			//The BC's change in velocity due to collisions
	double dL[3];			//The BC's change in angular velocity due to collisions

	// Qualities of the object
	int DSPLC;				// Flags whether bc can be pushed around by particles (0-no; 1-yes)
	double MASS;						//The BC's mass (only relevent if it moves)
	double VOL;				//Body's volume
	double I[3][3];		//The body's moment of inertia
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
typedef struct cell {
	int POP;					//Total population of the cell
	double MASS;					//Total mass of the cell
	int SP[MAXSPECI];				//Subpopulations of each species type in the cell
	double I[3][3];		//The cell's moment of inertia about (0,0,0)
	double E[3][3];		//Velocity gradient tensor --- first index [i] is on velocity, second [j] on derivative --- E[i][j]= dv[i]/dx[j]
	double S;					//The cell's scalar order parameter
	double CM[3];			//Centre of mass
	double VCM[3];		//Centre of mass velocity of the cell
	double DIR[3];		//Director of the cell (average orientation)
	double FLOW[3];		//Centre of mass velocity of the cell averaged over FLOWOUT time steps
	double Ps[3][3];		//Streaming part of the local instantaneous stress tensor
	double Pc[3][3];		//Collisional part of the local instantaneous stress tensor

	struct particleMPC *pp;		//Pointer to first SRD particle in list
	struct particleMD *MDpp;	//Pointer to first MD particle in list
	struct smono *sp;					//Pointer to first swimmer monomer particle in list
} cell;
typedef struct outputFilesList {
	FILE *fcoarse,*fflow,*fvel,*fenergy,*fenergyfield,*fenneighbours;
	FILE *fsynopsis,*favvel,*forder,*forderQ,*forderQK,*favs,*fdensSTD,*fchckpnt,*fenstrophy,*fmultiphase,*fpressure;
	FILE *fcorrVV,*fcorrNN,*fcorrWW,*fcorrDD,*fcorrSS,*fcorrPP,*fbinder;
	FILE *fhistVel,*fhistSpeed,*fhistVort,*fhistEnstr,*fhistDir,*fhistS,*fhistDens;
	FILE *fenergyspect,*fenstrophyspect;
	FILE *fdefects;
	FILE *fdetail[MAXSPECI];
	FILE *fsolids[MAXBC];
	FILE *fswimmers,*fswimmersOri,*fruntumble;
} outputFilesList;
typedef struct outputFlagsList {
	int TRAJOUT;				//Flag for if the detailed trajectories of every particle are outputted
	int COAROUT;				//Flag for if coarse grain is outputted
	int AVVELOUT;				//Flag for if total average velocity is outputted
	int FLOWOUT;				//Flag for if the flow field is outputted
    int VELOUT;                 //Flag for if the velocity field is outputted
	int HISTVELOUT,HISTSPEEDOUT,HISTVORTOUT,HISTENSTROUT,HISTDIROUT,HISTSOUT,HISTNOUT;	//Flag for if distributions are outputted
	int ENERGYSPECTOUT,ENSTROPHYSPECTOUT;	//Flag for if energy and enstrophy spectra are outputted
	int DEFECTOUT;			//Flag for if defect positions are outputted
	int ENOUT;					//Flag for if system energy is outputted
	int ENFIELDOUT,ENNEIGHBOURS;	//Flag for if orientational energy as a function of position is outputted
	int SPOUT;					//Flag for if the colour/phi/species-type field is outputted
	int PRESOUT;					//Flag for if the pressure field is outputted
	int SYNOUT;					//Flag for if system synopsis is outputted
	int ORDEROUT;				//Flag for if the order parameter are outputted
	int QTENSOUT,QKOUT;		//Flag for if the order parameter tensor (and reciprocal space Q) is outputted
	int AVSOUT;					//Flag for if total average scalar order parameter is outputted
	int DENSOUT;				//Flag for density standard deviation is outputted
	int ENSTROPHYOUT;		//Flag for if total average enstrophy is outputted
	int CHCKPNT;				//Flag for checkpointing
	int CHCKPNTrcvr;		//Flag for simulation from recovery of checkpoint
	int BINDER,BINDERBIN;	//Flag for Binder cumulant and the bin size of the Binder cumulant
	int printSP;				//How many of the species are printed
	int SOLOUT;					//How often the solid bc trajectories are outputted
	int CVVOUT,CNNOUT,CWWOUT;	//Flag for if velocity-vel, director-dir and vorticity-vort spatial correlations are outputted
	int CDDOUT,CSSOUT,CPPOUT;	//Flag for if density-dens, order-order and phi-phi (phase) spatial correlations are outputted
	int SWOUT,SWORIOUT,RTOUT;				//Flag for if swimmer positions, orientations or run/tumbles are outputted
} outputFlagsList;
typedef struct kinTheory {
	double MFP;					//Mean Free Path
	double VISC;				//Kinematic viscosity
	double THERMD;			//Thermal Diffusion
	double SDIFF;				//Self diffusion
	double SPEEDOFSOUND;	//Speed of sound
	double sumM;				//Sum of all masses
} kinTheory;
typedef struct inputList {
	double KBT;					//Temperature: A third of thermal energy --- sets energy scale
	double dt;					//MPCD time step value
	int stepsMD;				//Number of MD steps per SRD step (NOTE make variable)
	int warmupSteps,simSteps;	//Number of iterations in the warmup and simulation phases
	unsigned long seed;	//seed for random number generator
	double RA,C,S;			//Rotation Angle,cos(RA),sin(RA)
	double GRAV[_3D];		//Constant acceleration from external force
	double MAG[_3D];		//Constant external magnetic field to torque nematogens
	int GRAV_FLAG,MAG_FLAG;		//Flags if no acceleration or torque
	double FRICCO;			//Friction coefficient for Langevin thermostat
	int TSTECH;					//Temperature scaling technique
	double TAU;					//The temperature relaxation time scale
	double MFPOT;				//The mean-field potential from self-consistent mean-field liquid crystals
	int RTECH;					//Rotation technique
	int LC;							//If LC=LCG=2 then liquid crystal using global S, if LC=LCL=1 then use local S, else isotropic (ISOF=0)
	int RFRAME;						//Flags initial galilean trans to rest frame (0 No shift, 1 shift)
	int zeroNetMom;			//If GAL=1 AND GRAV[D3] = [0,0,0] then zero every zeroNetMom steps.
	int GALINV;					//If rshift=1 do the random shifting of mpcd particles. If zero then do not zeroNetMom steps.
	int MFPLAYERH; 				// Height above which MFP goes to 0. For simulating thin ordered films below disordered fluids. If 0, disable this functionality.
} inputList;
typedef struct specSwimmer {
	int TYPE;						//Type of swimmer 0=fix-dipole; 1=dumbbell; 2=dumbell with excluded vol; 3=dumbell with not counter-force
	int QDIST;					//QDIST is the flag which tells how the swimmers are intially placed
	int ODIST;					//ODIST is the flag which tells how the swimmers are intially placed
	int headM,middM;					//Mass of head and middle monomers in swimmer (tail mirrors the head)
	double iheadM,imiddM;			//Inverse masses
	int HSPid,MSPid;				//Multiphase fluid particle type of the head and middle monomers in the swimmer  (tail is invisible)
	double FS;					//Magnitude of propolsion force to set swimming speed
	double TS;					//Magnitude of torque to set swimmer rotlet dipole (sign sets CW or CCW)
	double DS;					//Dipole strength (DS>0 --->pusher; DS<0 --->puller) (i.e. length of dipole in multiples of head/middle separation)
	double sizeShrink;	//How much LJ sigma and ro are shrunk when tumbling
	double springShrink;//How much the spring constant is shrunk when tumbling
	double fixDist;			//The fixed distance from the wall for DUMBBELL_NEARWALL mode
	double k;						//Spring constant
	double ro,iro;			//Spring separation
	double sig,isig;		//LJ sigma --- diameter approx sigma
	double eps;					//LJ interaction energy
	double runTime,tumbleTime;	//Average run and tumble times (if either is zero then no run/tumble dynamics). In units of MPC time steps (dt)
	int shrinkTime;			//Set time to shrink/extend in units of MPC time steps (dt)
	double MAGMOM;			//Magnetic moment strength/magnitude
} specSwimmer;
typedef struct smono {
	double Q[_3D];			//Q is the position of the monomer
	double V[_3D];			//V is the velocity
	double A[_3D];			//A is the acceleration
	int HorM;						//Signfies whether monomer is a head, or a middle (0=head; 1=middle) --- tail takes the mass of the head and doesn't participate
	struct smono *next;	//pointer to next particle in cell list
	struct smono *previous;	//pointer to previous particle in cell list
} smono;
typedef struct swimmer {
	smono H,M;					//The head, middle and tail of the dumbbell
	int RT;							//Flag for if in Run phase (0) or Tumble phase (1) or shrink phase (2) or extension phase (3)
	double n0[_3D];			//Orientation vector of the swimmer at start of phase
	int timeCNT;				//Time counter for run or tumble
	int timeRND;				//The randomly generated time that a given run or tumble will go for
	double ro,iro;			//FENE separation --- Each swimmer needs it's own since it changes for run/tumbler
	double sig,isig;		//LJ sigma --- diameter approx sigma --- Each swimmer needs it's own since it changes for run/tumbler
	double k;			//FENE separation --- Each swimmer needs it's own since it changes for run/tumbler
} swimmer;

#endif
