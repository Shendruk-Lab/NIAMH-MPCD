
//================================================================================
//
// name:   mdtypes.h
// author: ftessier
// date:   2005-05-02 @ 16:33:34
//
//================================================================================


#ifndef MDTYPES_H
#define MDTYPES_H


#include <stdio.h>
#include <sys/types.h>

#include "../mpcd/headers/SRDclss.h"

//================================================================================
// Macros
//================================================================================

// define these in the Makefile instead, for convenience
// #define PHASE_COUNT			2					// number of simulation phases
// #define TYPE_COUNT			3					// number of particle types
// #define POLY_SETS			5					// number of polymer sets
// #define Q_SETS				4					// number of polymer sets
// #define THERMOSTAT_DPD
// #define TEMPERATURE_CONF

// prototype for snprintf on MACI
#ifdef MACI
extern int snprintf (char *__restrict __s, size_t __maxlen, __const char *__restrict __format, ...);
#endif

// step counters for simulation phases
#define	PC					PHASE_COUNT+1		// shorthand notation for PHASE_COUNT+1
#define TC					TYPE_COUNT			// shorthand notation for TYPE_COUNT
#define	PS					POLY_SETS			// shorthand notation for POLY_SETS
#define	QS					Q_SETS				// shorthand notation for Q_SETS
#define count_					0

// version control
#define AUTHOR					"Frederic Tessier"

// parameter types
#define INTG   					0
#define HEXA   					1
#define REAL  					2
#define real					float
#define	real_FORMAT_STR 			"%f"

// coordinates
#define x_					0
#define y_					1
#define	z_					2
#define	xy_					3
#define	xz_					4
#define	yz_					5
#define xyz_					6
#define r_					0
#define q_					1
#define rq_					3
#define rz_					4
#define qz_					5
#define rqz_					6
#define p_					2
#define	rp_					4
#define	qp_					5
#define	rqp_					6

// coordinate systems
#define CARTESIAN				0x00100
#define	CYLINDRICAL				0x01000
#define	SPHERICAL				0x10000

//MD-MPC coupling modes
# define noMD 					0
# define MDinMPC 				1
# define MPCinMD 				2

// domains
#define	DOMAIN_ALL				-1

// measurement actions
#define ACTION_RESET				0
#define	ACTION_COLLECT				1
#define ACTION_AVERAGE				2
#define ACTION_NORMALIZE			3
#define	ACTION_PRINT				4
#define	ACTION_COMPUTE				5
#define	LOOP_ATOM				100
#define	LOOP_NEBR				101
#define	LOOP_POLYMER				102

// normalization actions
#define	NORM_NONE				-1
#define	NORM_AVERAGE				0
#define	NORM_MEASURE				1
#define	NORM_COUNT				2
#define	NORM_USER				3
#define	NORM_INTEGRAL				4
#define	NORM_VOLUME				5

// crystal lattice
#define LATT_SC					0
#define	LATT_BCC				1
#define	LATT_FCC				2
#define LATT_BULK				3

// predefined LJ cutoff
#define	CUTOFF_WCA				-1

// system geometry
#define GEOM_BULK				0
#define	GEOM_CAPILLARY				1
#define	GEOM_PLATES				2
#define	GEOM_CYLINDER				3

// density types
#define	VOLUME					0
#define	SURFACE					1

// layout rules
#define	LAYOUT_NONE				-1
#define	LAYOUT_WALL				0
#define	LAYOUT_FLUID				1
#define	LAYOUT_SURFACE			2
#define	LAYOUT_ANCHOR				3
#define	LAYOUT_PLATES				4
#define	LAYOUT_CYLINDER			5
// Tyler added the following
#define	LAYOUT_RODX				6
// Zahra added the following
#define	LAYOUT_RODY				7
// Karolina added the following
#define	LAYOUT_U				8

// atom types (index)
#define TYPE_WALL				0
#define	TYPE_FLUID     		1
#define TYPE_MONOMER				2

// atom groups (bitfield)
#define GROUP_MAX				30					// maximum number of groups
#define GROUP_NONE				0x0
#define GROUP_FLUID				0x00000001
#define GROUP_WALL				0x00000002
#define GROUP_WALL_IN				0x00000004
#define GROUP_WALL_OUT				0x00000008
#define GROUP_ION				0x00000010
#define GROUP_ION_POS				0x00000020
#define GROUP_ION_NEG				0x00000040
#define GROUP_MONOMER				0x00000100
#define GROUP_GRAFT				0x00000200
#define GROUP_ALL      				0x0FFFFFFF

// atom visual properties
#define	OBJ_FLUID				0
#define	OBJ_WALL				1
#define	OBJ_ION_POS				2
#define	OBJ_ION_NEG				3
#define	OBJ_PBCBOX				4
#define	OBJ_MONOMER				5
#define	OBJ_MONOMER_1				6
#define	OBJ_MONOMER_GRAFT 			7
#define	OBJ_MONOMER_GRAFT_1 			8

// pbc offset indices
#define PBC_INDEX_NONE 				13
#define PBC_INDEX_x 				9
#define PBC_INDEX_y 				3
#define PBC_INDEX_z 				1

// pbc condition bitfield
#define PBC_COND_x				0x001
#define	PBC_COND_y				0x002
#define	PBC_COND_z				0x004

// simulation options
#define	SIMOPT_SETUP_NEW			0
#define	SIMOPT_SETUP_CHKPOINT			1

// atom picking
#define	PICK_POINTER				0
#define PICK_RANDOM				-1
#define	PICK_LAST				-2
#define	CHECK					1
#define	noCHECK					0

// loop limits
#define	LOOP_MAX				100000				// to prevent infinite loops
#define	GROWLOOP_MAX				100				// to control polymer growth recursion

// number of dimensions (must be 3)
#define DIM_MD						3			// just to clarify the code

// varia
#define pi					M_PI
#define SQRT_3					1.732050807568877
#define	STRLEN					256
#define STRMAX					255
#define DONE					1
#define	YES					1
#define	NO					0
#define	NONE					-1


//================================================================================
// Macro functions
//================================================================================

// error and log messages
#define error(e)  			ReportError (e, __FILE__, __LINE__);
#define LOG(...)    			{ snprintf (sim->msg, STRMAX, __VA_ARGS__); Report (sim, GetSimStream(sim->files,"log"), sim->msg); }

// parameter initialization trick (from Rapaport)
#define  INTG_PARAM(x,y)		{ #y, &(x->y), INTG,  sizeof(x->y)/sizeof(int), 			0}
#define  REAL_PARAM(x,y)		{ #y, &(x->y), REAL,  sizeof(x->y)/sizeof(real), 			0}
#define  HEXA_PARAM(x,y)		{ #y, &(x->y), HEXA,  sizeof(x->y)/sizeof(unsigned int), 	0}


//================================================================================
// Type definitions
//================================================================================

// Try to fit the particle in a multiple of 32 bytes (the size of the L1 cache
// line on the SUN machines. The r, v, a, and t triplets are also split to each
// be on independent cache lines (yielding a speed gain of about 10%!). Also
// keep acceleration and mass on same cache line. The prev and next pointers are
// used exclusively for fene bonds in linear polymers. In time, this should be
// replaced by a general bond list, with bond type etc.


// simulation particle structure
typedef struct particleMD {				// a particle					double		 float
    real  	rx, ry, rz;					// position							24			12
    real  	kspring;					// spring constant when anchored 	 8		 	 4
    real  	vx, vy, vz;					// velocity							24			12
    int		anchor;						// anchor flag						 4			 4
    real	q;							// charge							 4			 4
    real  	ax, ay, az;					// acceleration						24			12
    real  	mass;						// mass								 8			 4
    real  	x0,  y0, z0;				// anchor position					24			12
    int    	type;						// type								 4			 4
    int		group;						// group bit field					 4		 	 4
    real  	wx, wy, wz;					// real world position				24			12
    real  	Tfx, Tfy, Tfz, Tdivf;		// for configurational temp.		32			16
    int		pad[7];						// memory padding					28			28
    int		object;						// object id (for viewer)			 4			 4
    struct  particleMD *prev, *next;		// previous and next monomer		 8			 8
    struct  particleMD *prevSRD, *nextSRD;	// previous and next monomer for SRD binning	 8			 8
} particleMD;	 							// TOTAL						   240


// point in 3d space
typedef struct point {
    int		coord;						///< coordinate system
    real  	c1, c2, c3;					///< coordinate in 3 dimensions
    real  	x0, y0, z0;					///< coordinate system cartesian origin
} point;


// simulation parameter structure
typedef struct simparam {
    char	*name;						///< name
    void	*value;						///< pointer to storage space
    int		type;						///< parameter type
    int		count;						///< vector length
    int		n;							///< number of values read
} simparam;
typedef simparam *paramptr;


// simulation file structure
typedef struct simfile {
    FILE	*stream;					///< file descriptor
    char	label[STRLEN];				///< file
    char	name[STRLEN];				///< file name
    char	desc[STRLEN];				///< file description
    char	type[STRLEN];				///< content type (ascii, binary, ...)
    char	layout[STRLEN];				///< content layout
    char	columns[STRLEN];			///< column names
    int		frame;						///< current frame index
    long	chkpos;						///< file position of last checkpoint
    struct simfile *next;				///< pointer to the next file in list
} simfile;
typedef simfile *fileptr;				///< pointer to a file structure


// simulation scene structure
typedef struct simscene {

	// counters
	int			stepPrint[PC];			///< interval between scene prints
	int			frame;					///< current frame index

	// scene function pointers
	void		(*scenefunc) (void *sim, struct simscene *scene, int action);

	// scene information
	char		label[STRLEN];			///< scene label
	int			active;					///< scene is active?
	int			file;					///< scene linked to a file?
	FILE		*stream;				///< associated data stream

	// atom group filters
    int			groupInc;				///< included atom group(s)
    int			groupExc;				///< included atom group(s)

    // space domain
    int			sd;						///< space domain dimensions
    int			scoord;					///< space domain coordinate system
    int			sdomain;				///< space domain
    real  		s1min, s2min, s3min;	///< domain minima
    real  		s1max, s2max, s3max;	///< domain maxima

	// scene pointers
    struct simscene *next;				///< pointer to the next scene
} simscene;
typedef simscene *sceneptr;				///< pointer to a scene structure


// simulation histogram strucutre
typedef struct simhist {

    // counters
    int			stepCollect[PC];		///< interval between data collection
    int			stepPrint[PC];			///< interval between prints
    int			frame;					///< current frame index

    // histogram function pointers
    void    	(*histfunc)  (void *sim, struct simhist *h, int action);
    int			(*func)      (void *sim, struct simhist *h, int action);

    // histogram function names
    char		histfuncName[STRLEN];	///< histogram function name
    char		funcName[STRLEN];		///< histogram value function name

    // histogram information
    char		label[STRLEN];			///< histogram label
    int			active;					///< data collection (on/off)
    int			file;					///< histogram linked to file (yes/no)
    FILE		*stream;				///< associated data stream

    // atom group filters
    int			groupInc;				///< included atom group(s)
    int			groupExc;				///< included atom group(s)

    // histogram data
    real  		*bin, *count;			///< data bins and bin element counts
    real  		totCount;				///< total element count
    int			n;						///< total number of bins
    int			m;						///< number of measurements
    real  		norm;					///< normalization factor
    int			normMode;				///< which type of normalization to use

    // measurement variables
    int			histLoop;				///< type of loop
    particleMD	*p1, *p2;				///< pointers to particles
    point		*r, *rh;				///< pointer to points
	real		*hc1, *hc2, *hc3;		///< pointer to coordinates
	real		value;					///< value to add in histogram

    // histogram domain
    int			d;						///< histogram dimension (1,2,3)
    int			coord;					///< coordinate system
    int			domain;					///< histogram domain
    int			n1,   n2,   n3;			///< number of bins
    real  		min1, min2, min3;		///< minima
    real  		max1, max2, max3;		///< maxima
    real  		bin1, bin2, bin3;		///< bin sizes

    // space domain
    int			sd;						///< space domain dimensions
    int			scoord;					///< space domain coordinate system
    int			sdomain;				///< space domain
    real  		smin1, smin2, smin3;	///< domain minima
    real  		smax1, smax2, smax3;	///< domain maxima

    // histogram pointers
	struct simhist *h1;					///< friend histogram
    struct simhist *next;				///< pointer to the next hist

} simhist;
typedef simhist *histptr;				///< pointer to a hist structure


// simulation options
typedef struct simoptions {
    int setupType;						///< type of simulation setup
} simoptions;


// elements of the list1STD list
typedef struct item1STD {
    particleMD	*p1;					///< pointer to a particle
} item1STD;


// elements of the list2STD list
typedef struct item2STD {
    particleMD	*p1, *p2;				///< pointers to the two particles
} item2STD;


// elements of the list3STD list
typedef struct item3STD {
    particleMD	*p1, *p2, *p3;				///< pointers to the three particles
} item3STD;


// elements of the list2PBC list
typedef struct item2PBC {
    particleMD	*p1, *p2;				///< pointers to the two particles
    real  		*pbc;					///< pre-computed PBC offset for this pair
} item2PBC;


// elements of the listPoly list
typedef struct	itemPoly {
	particleMD	*p1;					///< pointer to the first monomer
	int			grafted;				///< grafted flag
} itemPoly;


// elements of the listCount list
typedef struct	itemStep {
    int			*counter;				///< pointer to actual counter variable
    void		(*action) (void *sim);	///< action function pointer
} itemStep;


// list of single particle pointers
typedef struct list1STD {
	item1STD	*items;					///< items of the list
	int			n, max;					///< current and maximum number of items
} list1STD;


// list of pairs of std particles
typedef struct list2STD {
  item2STD	*items;					///< items of the list
	int			n;						///< current number of items in the list
	int			max;					///< maximum number of items
} list2STD;

// list of triplets of std particles
typedef struct list3STD {
  item3STD	*items;					///< items of the list
	int			n;						///< current number of items in the list
	int			max;					///< maximum number of items
} list3STD;

// list of pairs of pbc particles
typedef struct list2PBC {
	item2PBC	*items;					///< items of the list
	int			n;						///< current number of items in the list
	int			max;					///< maximum number of items
} list2PBC;


// list of polymers
typedef struct listPoly {
	itemPoly	*items;					///< items of the list (the polymers here)
	int			n;						///< current number of items in the list
	int			max;					///< maximum number of items
} listPoly;


// list of atoms
typedef struct listAtom {
	particleMD	*items;					///< items of the list (particles here)
	int			n;						///< current number of items in the list
	int			max;					///< maximum number of items
} listAtom;


// list of counters
typedef struct listStep {
	itemStep	*items;					///< items of the list
	int			n;						///< current number of items in the list
	int			max;					///< maximum number of items
} listStep;


//================================================================================
// Main simulation data structure type definition
//================================================================================
typedef struct simulation {		 		// a simulation

    // atom lists
    listAtom	atom;			 		///< list of all atoms
    list1STD	anchor;		 			///< list of anchored atoms
    list1STD	charge;		 			///< list of coulombic charges
    list2STD	fene;					///< list of fene interaction pairs
    list3STD	bend;					///< list of bend interaction triplets
    listPoly	polymer;				///< list of polymers

    // neighbors
    list2STD	nebrSTD;		 		///< neighbor pair pointers list (STD for standard)
    list2PBC    nebrPBC;		 		///< neighbor pair pointers list (PBC for periodic boundary conditions)
    real  		rNebrShell;		 		///< width of the neighbor region shell
    real  		rNebr;			 		///< neighborhood radius (rCut + rNebrShell)
    real  		rNebr2;			 		///< neighborhood radius squared (rNebr^2)
    real  		**pbc;			 		///< periodic boundary offsets
    int			pbcond;			 		///< periodic boundary condition
    real  		drTotMax;				///< maximum displacement (controls nebr list refresh)

    // neighbor cells
    int 		***cell;		 		///< list of cell entrypoints
    int			*cellList;		 		///< list of atoms linked into cells
    int			nCell;			 		///< total number of cells
    int			nCellAxis[DIM_MD];		 	///< number of cells along x, y, z
    real  		cellInvWidth[DIM_MD];	 	///< 1 / (width of cells) along x, y, z

    // polymers
	int			polyAtomType[PS];		///< type for polymer atoms
	int			polyLayout[PS];			///< where to put the polymers (e.g., BULK, SURFACE, ...)
	int			polyDensityKind[PS]; 	///< kind of density provided (e.g., BULK, SURFACE, ...)
	real		polyDensity[PS];		///< density of polymers
	real		polySpread[PS];			///< spacing factor between SURFACE polymers
    int			polyM[PS];				///< number of polymer chains
    int			polyN[PS];				///< number of monomers per chain
	real			monoCharge[PS];			///< charge of each monomer

	// charges
	int			qLayout[QS];  			///< where to put the charges (e.g., SURFACE, TYPE_FLUID, ...)
	int			qDensityKind[QS];		///< kind of density provided (e.g., BULK, SURFACE, ...)
	real		qDensity[QS];			///< the density of charges
	real		qSpread[QS];			///< spacing factor for SURFACE charges
	int			qCharge[QS];			///< the charge in units of the elementary charge
	int			qNumber[QS];			///< the number of charges to put (will be compounded with density)

	// layout quantities
    int			nPore;					///< number of sites in capillary pore
    real   		poreVolume;				///< pore volume from void wall sites
	real		surfaceArea;			///< surface area
	list1STD	wall;					///< list of wall atoms
	list1STD	fluid;					///< list of fluid atoms
	list1STD	surface;				///< list of surface atoms

    // time step
    real  		tNow;			 		///< current md time
    real      	dt;						///< md timestep increment

    // relaxation parameters
    int			stepRelax;				///< relaxation flag
    real  		drRelax;				///< relaxation position increment

    // physical values
    real		kT[PC];					///< thermal energy in units of epsilon
    real  		rCut;			 		///< cutoff distance of Lennard-Jones interaction
    real  		ljShift;		 		///< LJ energy shift
    real		sigma_lj;
    real		lambda_D;
    real  		rCutCoul;		 		///< coulomb interaction cutoff
    real  		bjerrum;		 		///< bjerrum length of the fluid
    real  		condenseCriteria;	 		///< condensation criteria for mpc particles
    real  		r0Fene;		 			///< maximal extension of FENE bonds
    real  		kFene;		 			///< strength of the FENE potential
    real  		kSqu;		 			///< strength of the 2D harmonic potential
    real  		kBend;		 			///< strength of the bend potential
    real  		kNemMPC;		 			///< strength of the bend potential with the background nematic MPCD
    real  		Efield[PC];				///< external electric field
    real  		overlapMin;				///< min initial distance between atoms
    real  		overlapMinMonomer;		///< min initial distance between monomers

    // atom type parameters
    real		rhoType[TC];			///< dimensionless density of particles
    int			nType[TC];	 			///< number of atoms
    int			typeStart[TC];	 		///< first atom
    real		atomMass[TC];	 		///< mass
    int			atomAnchor[TC];	 		///< anchor flag
    real  		atomSpring[TC];	 		///< spring constant (if anchored)
    real   		qPosNum[TC];			///< number of positive ions
    real  		qNegNum[TC];			///< number of negative ions
    real  		qPosRho[TC];			///< volume density of poitive ions
    real  		qNegRho[TC];			///< volume density of negative ions
    real  		qPosSigma[TC];			///< surface density of negative ions
    real  		qNegSigma[TC];			///< surface density of negative ions
    int			qPosValence[TC];		///< valence of positive ions
    int			qNegValence[TC];		///< valence of negative ions
    real  		qLayer[TC];	 			///< thickness of charged layer

    // step counters
	listStep	stepCounter;			///< list of step counters
    int			phase;			 		///< current simulation phase
    int			nStep[PC];				///< number of steps in each phase
    int			step[PC];				///< accumulated number of steps
    int			stepAvg[PC];			///< steps between averages
    int			stepCheckpoint[PC];		///< steps between checkpoints


    // histogram step counters
    int			stepHistCollect[PC];	///< interval between histogram collection
    int			stepHistPrint[PC];		///< interval between histogram writes

	// scene step counters
    int			stepSceneSlow[PC];		///< steps between scene output
    int			stepSceneFast[PC];		///< steps between scene output

    // measurements
    real  		totE;			 		///< total energy
    real  		kinE;			 		///< kinetic energy
    real  		potE;			 		///< total potential energy
    real  		ljE, harmE;				///< potential energies
    real  		coulE, feneE;			///< potential energies
    real  		bendE, nemE;				///< potential energies
    real  		s_kinE, ss_kinE;	 	///< kinetic energy accumulators
    real  		s_potE, ss_potE;	 	///< potential energy accumulators
    real  		s_totE, ss_totE;	 	///< total energy accumulators
    real  		kinETherm;				///< kinetic energy for thermostat
    histptr		histograms;		 		///< histogram list
	sceneptr	scenes;					///< scene list

    // thermostat
    int			groupThermRescale[PC];	///< group of atoms subject to rescaling
    int			stepThermRescale[PC];	///< steps between temperature rescalings
    int			groupThermDPD[PC];		///< group of atoms subject to DPD thermostat
    real  		eta[PC];				///< DPD thermostat dissipation strength
    int			nAtomThermRescale;		///< number of atoms subject to rescaling
    int			nAtomThermDPD;			///< number of atoms subject to DPD

    // simulation box
    real       	unitCells[DIM_MD];		 	///< number of unit cells along each axis
    int			lattice;		 		///< what type of crystal lattice
    int			geometry;		 		///< system geometry
    real		nAtomCell;		 		///< number of atoms per unit cell
    real  		box[DIM_MD];		 		///< dimensions of the simulation box
    real  		boxHalf[DIM_MD];		 	///< 0.5*box
    real  		rho;					///< density of the simulation box

    // capillary
    real  		wallThickness;		 	///< wall thickness (in unit cells)
    real  		caprIn, caprIn2;		///< inside radius
    real  		caprOut, caprOut2;		///< ouside radius
    real  		caprInTrue;				///< unambiguous internal capillary radius
    real  		dLattice;		 		///< distance between n.n. lattice atoms

    // simulation parameters
    paramptr	param;		 			///< pointer to a list of all parameters
    int			nParam;					///< number of simulation parameters
	int			randomSeed;				///< Random number generator seed
    int			warmupMD;				///< = 0 No, = 1 yes
    // program information
    pid_t		pid;			 		///< process id of the simulation
    real		version;		 		///< version of the program

    // simulation files etc.
    simfile		*files;			 		///< linked list of data file pointers
    FILE		*input;					///< input file stream
    char		inputFile[STRLEN];		///< name of the input file
    char		outputDir[STRLEN];		///< directory where output files are written
    char		programName[STRLEN];	///< litteral name of the program
    char		label[STRLEN];			///< simulation label
    char		msg[STRLEN];		 	///< message buffer

} simulation;
typedef simulation *simptr;


#endif
