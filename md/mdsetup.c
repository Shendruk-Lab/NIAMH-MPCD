
//================================================================================
//
// name:   mdsetup.c
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
// Setup functions for the molecular dynamics program.
//
//================================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "mdtypes.h"
#include "mderror.h"
#include "mdutil.h"
#include "mdmemory.h"
#include "mdfiles.h"
#include "mdhistogram.h"
// #include "mdscene.h"
#include "mdvmd.h"
#include "mdmeasure.h"
#include "mdthermostat.h"
#include "mdsrd.h"
#include "mdsetup.h"

#ifdef MPI
#include <mpi.h>
#endif

/// Entry point for the initialization of the whole MD simulation. This
/// function parses the command-line arguments, allocates memory for the simulation
/// structure, and performs all initialization tasks. It initializes some fields of
/// sim and calls other functions to setup the parameters and the simulation world.
/// Contains conditional MPI compile branches to perform embarrassingly parallel
/// runs on multi-CPU systems. There are two main cases to handle: new simulations
/// and simulation started from a checkpointed state. Two separate functions take
/// care of these two cases.
///
/// @param		argc the number of command-line arguments
/// @param		argv the array of literal command-line arguments
/// @return		a pointer the newly created and initialized sim structure
/// @warning	At the end of this function, the simulation must be completely
///				ready to start, and the program immediately enters the main
///             MD step loop. Thus, all initialization must be included or
///				called from this function.

//================================================================================
simptr SetupSimulation (int argc, char *argv[])
//================================================================================
{
	simptr 		sim;
	simoptions	options;

	// MPI
	#ifdef MPI
	int  		i, nproc, rank;
	char 		chkname[128], *ptr;
	FILE 		*chklist;
	MPI_Status  status;
	#endif

	// allocate memory for the simulation information
	sim = (simptr) mycalloc (1, sizeof(simulation));

	// default options
	options.setupType = SIMOPT_SETUP_NEW;

	// parse command-line options
	ParseOptions (argc, argv, sim, &options);


	// MPI
	#ifdef MPI
	MPI_Init 		(&argc, &argv);
	MPI_Comm_size 	(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank 	(MPI_COMM_WORLD, &rank);

	#endif

	// ---------------------------------------------------------------------------
	// NEW simulation setup
	// ---------------------------------------------------------------------------
	if (options.setupType == SIMOPT_SETUP_NEW) {

		// MPI
		#ifdef MPI
		snprintf (sim->inputFile, STRMAX, "%s-%02d", sim->inputFile, rank);
		printf   ("MPI rank %d: %s\n", rank, sim->inputFile);
		#endif

		// setup the simulation parameters
		SetupParameters (sim);


		// setup simulation world
		SetupNewWorld (sim);

		// report
		LOG ("START of main simulation loop\n");
	}


	// ---------------------------------------------------------------------------
	// CHECKPOINTED simulation setup
	// ---------------------------------------------------------------------------
	else if (options.setupType == SIMOPT_SETUP_CHKPOINT) {

		// MPI
		#ifdef MPI
		// dispatch chk file names
		if (rank == 0) {
			// open chk file list
			if ((chklist = fopen (sim->outputDir, "r")) == NULL) error(EFILE);
			// rank 0 handles first file
			fgets (sim->outputDir, STRLEN, chklist);
			// rank i handle remaining files
			i=1;
			while (!feof(chklist)) {
				fgets (chkname, 128, chklist);
				MPI_Send (chkname, strlen(chkname)+1, MPI_CHAR, i, 99, MPI_COMM_WORLD);
				i++;
			}
			// close file list
			fclose (chklist);
		}
		else {
			MPI_Recv (sim->outputDir, STRLEN, MPI_CHAR, 0, 99, MPI_COMM_WORLD, &status);
		}
		// chomp and print name
		if (ptr = strchr(sim->outputDir,'\n')) *ptr = '\0';
		printf ("%d: %s\n", rank, sim->outputDir);
		#endif

		// read checkpoint file
		CheckpointRead (sim);

		// setup world from checkpoint
		SetupCheckpointWorld (sim);

		// report
		LOG ("START of main simulation loop\n");
	}

	// ---------------------------------------------------------------------------
	// UNDEFINED simulation setup
	// ---------------------------------------------------------------------------
	else {
		error (ESETUP);
	}

	// return a pointer to the simulation structure
	return sim;
}


/// Reads the simulation parameters. Required parameters are indicated in
/// the local array paramList. Calls the function ReadParameters() to handle
/// the actual scanning of the input.
///
/// @param		sim a pointer to a simulation structure
/// @return		void
/// @warning	Parameters added here must also be added to the sim
///				structure (and in the inpu file, obviously).

//================================================================================
void SetupParameters (simptr sim)
//================================================================================
{
	int   		n;
	simparam 	param[] = {

				// physical characteristics
				INTG_PARAM (sim, randomSeed),
				INTG_PARAM (sim, warmupMD),
				REAL_PARAM (sim, unitCells),
				INTG_PARAM (sim, lattice),
				INTG_PARAM (sim, geometry),
				REAL_PARAM (sim, rCut),
				REAL_PARAM (sim, rNebrShell),
				REAL_PARAM (sim, sigma_lj),
				REAL_PARAM (sim, lambda_D),
				REAL_PARAM (sim, bjerrum),
				REAL_PARAM (sim, rCutCoul),
				REAL_PARAM (sim, condenseCriteria),
				REAL_PARAM (sim, wallThickness),
				REAL_PARAM (sim, r0Fene),
				REAL_PARAM (sim, kFene),
				REAL_PARAM (sim, kSqu),
				REAL_PARAM (sim, kBend),
				REAL_PARAM (sim, kNemMPC),
				REAL_PARAM (sim, overlapMin),
				REAL_PARAM (sim, overlapMinMonomer),
				REAL_PARAM (sim, dt),
				REAL_PARAM (sim, drRelax),
				INTG_PARAM (sim, stepRelax),

				// particle types
				REAL_PARAM (sim, rhoType),
				REAL_PARAM (sim, atomMass),
				INTG_PARAM (sim, atomAnchor),
				REAL_PARAM (sim, atomSpring),

				// polymers
				INTG_PARAM (sim, polyAtomType),
				INTG_PARAM (sim, polyLayout),
				INTG_PARAM (sim, polyDensityKind),
				REAL_PARAM (sim, polyDensity),
				REAL_PARAM (sim, polySpread),
				INTG_PARAM (sim, polyM),
				INTG_PARAM (sim, polyN),
				REAL_PARAM (sim, monoCharge),

				// charges
				INTG_PARAM (sim, qLayout),
				INTG_PARAM (sim, qDensityKind),
				REAL_PARAM (sim, qDensity),
				REAL_PARAM (sim, qSpread),
				INTG_PARAM (sim, qCharge),
				INTG_PARAM (sim, qNumber),

				// data collection
				INTG_PARAM (sim, nStep),
				INTG_PARAM (sim, stepAvg),
				INTG_PARAM (sim, stepCheckpoint),

				// histogram step counters
				INTG_PARAM (sim, stepHistCollect),
				INTG_PARAM (sim, stepHistPrint),

				// scene step counters
				INTG_PARAM (sim, stepSceneSlow),
				INTG_PARAM (sim, stepSceneFast),

				// external field
				REAL_PARAM (sim, Efield),

				// temperature
				REAL_PARAM (sim, kT),
				HEXA_PARAM (sim, groupThermRescale),
				INTG_PARAM (sim, stepThermRescale),
				HEXA_PARAM (sim, groupThermDPD),
				REAL_PARAM (sim, eta)
	};

	// compute the number of parameters
	n = sizeof(param)/sizeof(simparam);

	// read the n parameters
	ReadParameters (sim->inputFile, param, n, "#parameter:", "#end-header");

	// copy parameters to the simulation structure
	sim->nParam = n;
	sim->param  = (paramptr) mycalloc (n, sizeof(simparam));
	memcpy (sim->param, param, n*sizeof(simparam));
}


/// Initializes a new simulation environment, such as the random number
/// generator, the simulation box size, etc. Finally calls a number of other
/// functions to initialize particle positions, velocities, etc.
///
/// @param		sim a pointer to a simulation structure
/// @return		void
/// @warning	Don't use this function when starting a simulation from a
///				checkpoint file.

//================================================================================
void SetupNewWorld (simptr sim)
//================================================================================
{

	int		  i;
	real	  *box, *boxHalf;

	// local sim variables
	box	 = sim->box;
	boxHalf = sim->boxHalf;

	// set initial simulation phase
	sim->phase = 1;
	while (sim->nStep[sim->phase] <= 0) sim->phase++;

	// set label, working directory, and data files
	SetSimLabel	(sim->label, sim->pid);
	SetWorkingDir  (sim->outputDir, sim->label);
	SetupDataFiles (sim);

	// report
	LOG ("--------------------------------------------------------------\n");
	LOG ("Setting up simulation world\n");

	// initialize the random number generator
	SetupRandomGenerator (sim);

	// set the LJ cutoff and shift for WCA potential
	if (sim->rCut == CUTOFF_WCA) {
		sim->rCut = pow (2, 1.0/6.0);
		sim->ljShift = 1;
	}

	// LATTICE dependent settings
	switch (sim->lattice) {
		case LATT_FCC:
			sim->nAtomCell = 4.;
			break;
		case LATT_BULK:
			sim->nAtomCell = 1.;
			break;
		default:
			error (ELATTICE);
	}

	// GEOMETRY depedent settings
	switch (sim->geometry) {
		case GEOM_CAPILLARY:
			sim->pbcond = PBC_COND_x;
			sim->rho	= sim->rhoType[TYPE_WALL];
			break;
		case GEOM_PLATES:
			sim->pbcond = PBC_COND_x | PBC_COND_z;
			sim->rho	= sim->rhoType[TYPE_FLUID];
			break;
		case GEOM_CYLINDER:
			sim->pbcond = PBC_COND_x;
			sim->rho	= sim->rhoType[TYPE_FLUID];
			break;
		case GEOM_BULK:
			sim->pbcond = PBC_COND_x | PBC_COND_y | PBC_COND_z;
			sim->rho = sim->rhoType[TYPE_FLUID];
			break;
		default:
			error (EGEOMETRY);
	}

	// calculate the initial lattice spacing
	sim->dLattice = pow( (sim->nAtomCell/ (sim->rho?sim->rho:1.) ), (1.0/3.0) ) / sqrt(2.0);

	// set the box size from rho and the specified number of unit cells
	for (i=0; i<DIM_MD; i++) {
// 		box[i]	 = sim->unitCells[i] / pow (sim->rho/sim->nAtomCell, 1.0/3.0);
		box[i]	 = sim->unitCells[i] * pow( (sim->nAtomCell/ (sim->rho?sim->rho:1.) ), (1.0/3.0) );
		boxHalf[i] = 0.5 * box[i];
	}

	// set the bulk volume;
	sim->poreVolume = box[x_]*box[y_]*box[z_];

	// report
	LOG ("  %-20s = %g\n",   "rCut",		 sim->rCut);
	LOG ("  %-20s = %g\n",   "ljShiftE",	 sim->ljShift);
	LOG ("  %-20s = %g\n",   "nAtomCell",	 sim->nAtomCell);
	LOG ("  %-20s = %g\n",   "dLattice",	 sim->dLattice);
	LOG ("  %-20s = (%f, %f, %f)\n", "box",  box[x_], box[y_], box[z_]);
	LOG ("  %-20s = (%f, %f, %f)\n", "unit Cells",  sim->unitCells[x_], sim->unitCells[y_], sim->unitCells[z_]);
	LOG ("  %-20s = %#05X\n", "pbcond",	  	 sim->pbcond);

	// setup step counter lists
	SetupStepCounters (sim);

	// setup particles
	SetupParticles (sim);

	// setup particle lists
	SetupPolymerList  (sim);
	SetupChargeList   (sim);
// 	SetupNeighborList (sim);
	SetupAnchorList   (sim);
	SetupFeneList	  (sim);
	SetupBendList	  (sim);

	// setup groups
	SetupGroups (sim);

	// setup measurements
	SetupMeasurements (sim);

	// jiggle and relax particle positions
//	JiggleParticles (sim);
	RelaxParticles  (sim);

	// set simulation time
	sim->tNow = 0;

	// start simulation phase
	SimulationPhaseReset (sim);

	// print initial scenes
	VMDPrint (sim, sim->scenes);

	// report
	LOG ("Simulation is ready to start\n");
	LOG ("--------------------------------------------------------------\n");
}


/// Initializes the simulation environment in the case we are starting from a
/// checkpointed simulation state. Most of the simulation state is read from
/// the checkpoint file, so there is little to do, except rebuild lists which
/// are not saved in the checkpoint.
///
/// @param		sim a pointer to a simulation structure
/// @return		void
/// @warning	Don't use this function when setting up a new simulation.

//================================================================================
void SetupCheckpointWorld (simptr sim)
//================================================================================
{

	// reopen data files
	ReopenDataFiles (sim);

	// report
	LOG ("--------------------------------------------------------------\n");
	LOG ("Setting up simulation world from checkpoint file\n");

	// initialize the random number generator
	SetupRandomGenerator (sim);

	// setup list of step counters
	SetupStepCounters (sim);

	// setup particle lists
	SetupPolymerList  (sim);
	SetupChargeList   (sim);
// 	SetupNeighborList (sim);
	SetupAnchorList   (sim);
	SetupFeneList	  (sim);
	SetupBendList	  (sim);

	// setup histogram streams and function pointers
	HistogramSetupFiles (sim);
	HistogramSetupFuncs (sim);

	// report
	LOG ("Simulation is ready to resume\n");
	LOG ("--------------------------------------------------------------\n");
}


/// Initializes the random number generator. Pass 0 to Random seed in order
/// to use the varying time(0)*getpid() as the seed; this is set by the
/// simulation parameter randomSeed.
///
///	@param		sim a pointer to a simulation structure
///	@return		void
/// @warning	If the randomSeed in sim is not set to 0, then the seed is fixed.
///				A warning is printed to the log file when this is the case.

//================================================================================
void SetupRandomGenerator (simptr sim)
//================================================================================
{
	unsigned long seed = RandomSeed(sim->randomSeed);
	LOG ("  %-20s = %lu\n", "seed", seed);

	// warn if we are using a fixed seed
	if (sim->randomSeed != 0) LOG ("  %s\n", "*** WARNING *** the random number generator seed is FIXED!");
}


/// Initializes all particles. Actually calls other functions to allocate memory,
/// set up particle positions and velocities, and to build the necessary
/// particle lists.
///
/// @param		sim a pointer to a simulation structure
/// @return		void
///	@warning	It is important to call GrowParticleList early, before we start
/// 			adding particles, because the AtomInsert function may try insertions
/// 			before the list size is checked.
/// @warning	Don't forget to set up the surface list when building a system
///				with surfaces!
/// @warning	Set up the polymers before the charges, because the former are more
///				difficult to insert owing to the random growth and spacing constraints.

//================================================================================
void SetupParticles (simptr sim)
//================================================================================
{
	// allocate a first memory block for the particles
	ResetListAtom (&sim->atom);
	GrowListAtom  (&sim->atom);

	// initialize particle positions according to geometry
	switch (sim->geometry) {

		case GEOM_BULK:
			InitRandomCoord  	(sim, TYPE_FLUID, LAYOUT_FLUID);
			SetupLayoutLists	(sim);
			break;

		case GEOM_PLATES:
			sim->surfaceArea	= sim->box[x_]*sim->box[z_];
			sim->poreVolume		= sim->box[x_]*sim->box[y_]*sim->box[z_];
			InitRandomCoord  	(sim, TYPE_FLUID, LAYOUT_FLUID);
			SetupLayoutLists	(sim);
			break;

		case GEOM_CYLINDER:
			sim->caprIn   = (sim->box[y_]-2.0)/2.;
			sim->caprOut  = (sim->box[y_]-2.0)/2.;
			sim->caprIn2  = (sim->box[y_]-2.0)*(sim->box[y_]-2.0)/4.;
			sim->caprOut2 = (sim->box[y_]-2.0)*(sim->box[y_]-2.0)/4.;
			LOG ("# capillary radius = %E\n", sim->caprIn);
			sim->surfaceArea	= M_PI*sim->box[x_]*(sim->box[y_]-2.0);
			LOG ("# surface Area = %E\n", sim->surfaceArea);
			sim->poreVolume		= M_PI*sim->box[x_]*(sim->box[y_]-2.0)*(sim->box[z_]-2.0)/4.0;
			LOG ("# Volume = %E\n", sim->poreVolume);
			InitRandomCoord  	(sim, TYPE_FLUID, LAYOUT_FLUID);
			SetupLayoutLists	(sim);
			break;

		case GEOM_CAPILLARY:
			InitCapillaryRadius	(sim);
			InitLatticeCoord   	(sim, TYPE_WALL, LAYOUT_WALL);
			InitCapillaryPore	(sim);
			InitRandomCoord		(sim, TYPE_FLUID, LAYOUT_FLUID);
			SetupLayoutLists	(sim);
			break;

		default:
			error (EGEOMETRY);
			break;
	}

	// initialize polymers and charges
	InitPolymers (sim);
	InitCharges  (sim);

	// initialize particle velocities
	InitVelocities (sim);

	// report
	LOG ("Particle initialization complete\n");
	LOG ("  %-20s = %d\n", "nAtom", sim->atom.n);
}


/// Creates particles of a given type by inserting them at random locations
/// inside the system (subject to layout rules and overlap constraints). The
/// number of atoms to add is first determined from the geometry and the
/// type, and the AtomInsert function is called for each insertion.
///
/// @param		sim a pointer to a simulation structure
/// @param		type the type of particles to insert
/// @param		layout the desired layout rule
/// @return		void
/// @warning	Note that particle pointers are NOT updated (noCHECK).

//================================================================================
void InitRandomCoord (simptr sim, int type, int layout)
//================================================================================
{
	int	toInsert;			// number of atoms to insert
	real	*unitCells;			// simulation cells
	real	poreVolume;			// volume of capillary pore
	real	*rhoType;			// densities
// 	real	nAtomCell;
	char	str[STRLEN];			// a string

	// local sim variables
	unitCells  = sim->unitCells;
	rhoType    = sim->rhoType;
	poreVolume = sim->poreVolume;
// 	nAtomCell  = sim->nAtomCell;

	// report
	LOG ("Initializing random atom coordinates\n");

	// initialize number of particles to insert
	toInsert = 0;

	// calculate how many particles to insert, according to geometry
	switch (sim->geometry) {

		case GEOM_BULK:
			switch (type) {
				case TYPE_FLUID:
//					toInsert = nAtomCell * unitCells[x_]*unitCells[y_]*unitCells[z_];
					toInsert = rhoType[type] * unitCells[x_]*unitCells[y_]*unitCells[z_];
					break;
				default:
					error (ETYPE);
			}
			break;

		case GEOM_PLATES:
			switch (type) {
				case TYPE_FLUID:
//					toInsert = nAtomCell * unitCells[x_]*unitCells[y_]*unitCells[z_];
					toInsert = rhoType[type] * unitCells[x_]*unitCells[y_]*unitCells[z_];
					if (toInsert<(sim->polyDensity[0]*sim->polyN[0]*sim->surfaceArea+sim->polyN[0]*sim->polyM[0])){
					  toInsert = sim->polyDensity[0]*sim->polyN[0]*sim->surfaceArea/sim->polyN[0];
					  toInsert = sim->polyM[0]*sim->polyN[0]+sim->polyN[0] * toInsert;}

					break;
				default:
					error (ETYPE);
			}
			break;
		case GEOM_CYLINDER:
				switch (type) {
				case TYPE_FLUID:
					toInsert = poreVolume * rhoType[type];
					if (toInsert<(sim->polyDensity[0]*sim->polyN[0]*sim->surfaceArea+sim->polyN[0]*sim->polyM[0])){
					  toInsert = sim->polyDensity[0]*sim->polyN[0]*sim->surfaceArea/sim->polyN[0];
					  toInsert = sim->polyM[0]*sim->polyN[0]+sim->polyN[0] * toInsert;}
					break;
				default:
					error (ETYPE);
			}
			break;

		case GEOM_CAPILLARY:
				switch (type) {
				case TYPE_FLUID:
					toInsert = poreVolume * rhoType[type];
					break;
				default:
					error (ETYPE);
			}
			break;

		default:
			error (EGEOMETRY);
	}

	// insert atoms
	while (toInsert--) AtomInsert (sim, type, layout, 0, noCHECK, noCHECK);

	// report
	TypeToStr (type, str);
	LOG ("  %-20s = %d atoms\n", str, sim->nType[type]);
}


/// Creates particles of a given type by placing them on lattice sites inside
/// the system (subject to the layout rules). The function iterates over
/// lattice sites and call the AtomInsert function for each new atom.
///
/// @param		sim a pointer to a simulation structure
/// @param		type the type of particles to insert
/// @param		layout the desired layout rule
/// @return		void
/// @warning	Only the FCC lattice has been implemented yet.
/// @todo		Implement other crystal lattices.

//================================================================================
void InitLatticeCoord (simptr sim, int type, int layout)
//================================================================================
{
	int			i, k, nx, ny, nz;		// counters
	real 		unitCells[DIM_MD], gap[DIM_MD];	// lattice cells
	real		*box, *boxHalf, rho;			// box dimensions
	real		rCell[DIM_MD], r[DIM_MD];		// cell and atom positions
	particleMD	thisAtom;				// particles
// 	particleMD	*atom;					// particles
	char		str[STRLEN];				// a string

	// local sim variables
// 	atom	= sim->atom.items;
	box	 	= sim->box;
	boxHalf = sim->boxHalf;
	rho	 	= sim->rhoType[type];

	// reset unit cell size (gap)
	for (i=0; i<DIM_MD; i++) {
		unitCells[i] = sim->unitCells[i];
		gap[i]  = box[i] / unitCells[i];
	}

	// handle densities different from the box rho
	if (rho != sim->rho) {
		for (i=0; i<DIM_MD; i++) {
			gap[i] *= pow(sim->rho/rho, 1.0/3);
			unitCells[i] = box[i] / gap[i];
			gap[i] = box[i] / unitCells[i];
		}
	}

	// report
	LOG ("Initializing lattice coordinates\n");
	LOG ("  %-20s = (%f, %f, %f)\n", "unit cells", unitCells[x_], unitCells[y_], unitCells[z_]);
	LOG ("  %-20s = (%f, %f, %f)\n", "cell gap", gap[x_], gap[y_], gap[z_]);

	// place particles on the lattice
	switch (sim->lattice) {
		case LATT_FCC:
			for (nx=0; nx < unitCells[x_]; nx++) {
				rCell[x_] = (nx+0.25) * gap[x_] - boxHalf[x_];
				for (ny=0; ny < unitCells[y_]; ny++) {
					rCell[y_] = (ny+0.25) * gap[y_] - boxHalf[y_];
					for (nz=0; nz < unitCells[z_]; nz++) {
						rCell[z_] = (nz+0.25) * gap[z_] - boxHalf[z_];
						// place 4 FCC atoms in unit cell
						for (i=0; i<4; i++) {
							for (k=0; k<DIM_MD; k++) {
								if (i==3 || i==k) r[k] = rCell[k];
								else r[k] = rCell[k]+0.5*gap[k];
							}
							// insert the atom
							thisAtom.rx = r[x_];
							thisAtom.ry = r[y_];
							thisAtom.rz = r[z_];
							AtomInsert (sim, type, layout, &thisAtom, noCHECK, noCHECK);
						}
					}
				}
			}
			break;
		default:
			error (ELATTICE);
	}

	// report
	TypeToStr (type, str);
	LOG ("  %-20s = %d atoms\n", str, sim->nType[type]);
}


/// Calculates the inner and outer capillary radii from the capillary
/// thickness parameter (wallThickness). The capillary is a cylinder with its
/// axis along the x-axis. The smallest of y and z is selected to calculate
/// the capillary diameter (i.e. tha capillary always has a circular cross-section).
///
/// @param		sim a pointer to a simulation structure
/// @return		void

//================================================================================
void InitCapillaryRadius (simptr sim)
//================================================================================
{
// 	int	y, z;
	real	*unitCells;
	real	*box, *boxHalf, cellSize;
	real	caprOut, caprOut2, caprIn, caprIn2;

	// local sim variables
	unitCells = sim->unitCells;
	box       = sim->box;
	boxHalf   = sim->boxHalf;

	// report
	LOG ("Initializing capillary radii\n");

	// calculate the capillary radii (leave one empty cell outside for safety)
// 	y = y_;
// 	z = z_;
// 	if (box[z_] <= box[y_]) {
// 		y = z_;
// 		z = y_;
// 	}
	cellSize = box[y_] / unitCells[y_];
	caprOut  = boxHalf[y_] - cellSize;
	caprIn   = caprOut - (sim->wallThickness * cellSize);
	caprOut2 = caprOut*caprOut;
	caprIn2  = caprIn*caprIn;

	// check that radii are consistent
	if (caprIn < 0 || caprIn > caprOut) error (ECAPILLARY);

	// save radii values in sim structure
	sim->caprIn   = caprIn;
	sim->caprOut  = caprOut;
	sim->caprIn2  = caprIn2;
	sim->caprOut2 = caprOut2;

	// report
	LOG ("  %-20s = %f\n", "caprOut", sim->caprOut);
	LOG ("  %-20s = %f\n", "caprIn",  sim->caprIn);
}


/// Calculates the number of missing wall atoms making up the pore. Along
/// with the density of the wall, this defines unambiguous values for the
/// capillary volume and its internal radius.
///
/// @param		sim a pointer to a simulation structure
/// @return		void
/// @warning	This function relies on the value of sim->nPore, which
///				must hold the count of wall atoms rejected to form the pore.
///				Thus, this variable must be correctly updated when building
///				the capillary wall; best is to put it in the LayoutRule function.

//================================================================================
void InitCapillaryPore (simptr sim)
//================================================================================
{
	// report
	LOG ("Initializing pore\n");

	// pore volume
	sim->poreVolume = sim->nPore/sim->rhoType[TYPE_WALL];
	LOG ("  %-20s = %d\n", "nPore",   sim->nPore);
	LOG ("  %-20s = %f\n", "poreVolume", sim->poreVolume);

	// true capillary radius
	sim->caprInTrue = sqrt (sim->poreVolume / (pi*sim->box[x_]));
	LOG ("  %-20s = %f\n", "caprInTrue",  sim->caprInTrue);

	// surface area
	sim->surfaceArea = 2*pi*sim->caprInTrue * sim->box[x_];
}


/// Insert polymers in the simulation. Fluid atoms are removed to leave room
/// for the monomers and thus keep the fluid density approximately constant. We then insert
/// the monomers in the form of randomly grown linear chains. The input file can
/// specify different polymer "sets", each with their own properties (e.g. length, number,
/// grafted or not, etc.). The function is a little involved to handle all the
/// different cases, so the reader is referred to comments in the code for more
/// details.
///
/// @param		sim a pointer to a simulation structure
/// @return		void
/// @warning		This function relies on the surface properties (list, area, etc.)
///				for grafted polymers, so make sure these are set properly

//================================================================================
void InitPolymers (simptr sim)
//================================================================================
{
	int		i, j, keep, grown, toRemove, set, c, picked, loop, toAdd;
	int		layout, surfaceMetric;
	int		polyBulkTot, polySurfaceTot;
	int		polyAtomType[POLY_SETS], polyLayout[POLY_SETS];
	int		polyDensityKind[POLY_SETS];
	int		polyM[POLY_SETS], polyN[POLY_SETS];
	list1STD	charge;
	real		monoCharge[POLY_SETS];
	real		polyDensity[POLY_SETS], polySpread[POLY_SETS];
	real		A, V, d2, d2Min;
	list1STD	*layoutList, *fluid, *surface, candidates;
	listPoly	polymer;
	particleMD	*p1, *p2, *p3, *atom;
	real		dx, dy, dz;
	real 		theta;

	// local sim variables
	A = sim->surfaceArea;
	V = sim->poreVolume;
	surface = &sim->surface;
	fluid = &sim->fluid;
	p3 = (particleMD *) mycalloc (1, sizeof(particleMD));

	// reset polymer and candidate lists
	ResetListPoly (&polymer);
	ResetList1STD (&candidates);

	// geometry variables
	surfaceMetric = CARTESIAN;
	switch (sim->geometry) {
		case GEOM_CAPILLARY:
			surfaceMetric = CYLINDRICAL | x_;
			break;
		default:
			break;
	}
	layout = NONE;
	layoutList = NULL;
	d2Min = 0;

	// polymer set variables
	polyBulkTot    = 0;
	polySurfaceTot = 0;
	for (set=0; set<POLY_SETS; set++) {

		polyLayout[set] = sim->polyLayout[set];
		if (polyLayout[set] == NONE) continue;

		polyAtomType[set]	= sim->polyAtomType[set];
		polyDensityKind[set]	= sim->polyDensityKind[set];
		polyDensity[set]	= sim->polyDensity[set];
		polySpread[set]	= sim->polySpread[set];
		polyM[set]		= sim->polyM[set];
		polyN[set]		= sim->polyN[set];
		monoCharge[set]	= sim->monoCharge[set];

		if (polyDensity[set] > 0) {
			switch (polyDensityKind[set]) {
				case SURFACE:	polyM[set] += polyDensity[set] * A;
								break;
				case VOLUME:	polyM[set] += polyDensity[set] * V;
								break;
				default:		error (EDENSITY);
								break;
			}
		}

		switch (polyLayout[set]) {
			case LAYOUT_SURFACE:	polySurfaceTot+=polyM[set];
									break;
			case LAYOUT_PLATES:		polySurfaceTot+=polyM[set];
									break;
			case LAYOUT_FLUID:		polyBulkTot+=polyM[set];
									break;
			// Tyler added this to place the polymer as a rod in x-direction - almost identical to LAYOUT_FLUID
			case LAYOUT_RODX:		polyBulkTot+=polyM[set];
									break;
			// Zahra added this to place the polymer as a rod in y-direction - almost identical to LAYOUT_RODX
			case LAYOUT_RODY:		polyBulkTot+=polyM[set];
									break;
			// Added by Karolina to initialise polymer as a hairpin
			case LAYOUT_U:		polyBulkTot+=polyM[set];
									break;
			case LAYOUT_ANCHOR:		polyM[set] = 1;
									polyBulkTot+=polyM[set];
									break;
			case LAYOUT_CYLINDER:	polySurfaceTot+=polyM[set];
									break;
			default:			error (ELAYOUT);
									break;
		}
	}

	// build polymer sets
	for (set=0; set<POLY_SETS; set++) {

		// skip empty polymer sets
		if (polyM[set] <= 0 || polyLayout[set]==NONE) continue;

		// report
		LOG ("Setting up polymer set %d:\n", set);

		// remove fluid atoms
		toRemove = polyM[set] * polyN[set];
		while (toRemove--) {
			if (!AtomRemove (sim, TYPE_FLUID, PICK_RANDOM, 0)) break;
		}
		// setup some variables, according to the layout
		switch (polyLayout[set]) {
			case LAYOUT_SURFACE:
				d2Min  	   = pow (polySpread[set]*sqrt(A/polySurfaceTot), 2.0);
				layout 	   = LAYOUT_FLUID;
				layoutList = surface;
				break;
			case LAYOUT_FLUID:
				d2Min  	   = pow (polySpread[set]*pow(V/polyBulkTot,1/3.0), 2.0);
				layout 	   = LAYOUT_FLUID;
				layoutList = fluid;
				break;
			// Tyler added - identical to LAYOUT_FLUID
			case LAYOUT_RODX:
				d2Min  	   = pow (polySpread[set]*pow(V/polyBulkTot,1/3.0), 2.0);
				layout 	   = LAYOUT_FLUID;
				layoutList = fluid;
				break;
			// Zahra added - identical to LAYOUT_RODX 
			case LAYOUT_RODY:
				d2Min  	   = pow (polySpread[set]*pow(V/polyBulkTot,1/3.0), 2.0);
				layout 	   = LAYOUT_FLUID;
				layoutList = fluid;
				break;
			case LAYOUT_U:
				d2Min  	   = pow (polySpread[set]*pow(V/polyBulkTot,1/3.0), 2.0);
				layout 	   = LAYOUT_FLUID;
				layoutList = fluid;
				break;
			case LAYOUT_ANCHOR:
				layout 	   = LAYOUT_FLUID;
				layoutList = NULL;
				break;
			case LAYOUT_PLATES:
				d2Min  	   = pow (polySpread[set]*sqrt(A/polySurfaceTot), 2.0);
				layout 	   = LAYOUT_FLUID;
				layoutList = NULL;
				break;
			case LAYOUT_CYLINDER:
				d2Min  	   = pow (polySpread[set]*sqrt(A/polySurfaceTot), 2.0);
				layout 	   = LAYOUT_CYLINDER;
				layoutList = NULL;
				break;
			default:
				error (ELAYOUT)
				break;
		}

		// rebuild candidate sites list (skipping busy atoms)
		if (candidates.items) free (candidates.items);
		ResetList1STD(&candidates);

		switch (polyLayout[set]) {
			case LAYOUT_ANCHOR:
				p1 = (particleMD *) mycalloc (1, sizeof(particleMD));
				AddItem1STD (&candidates, p1);			// register candidate (just one!)
				break;
			case LAYOUT_PLATES:
				p1 = (particleMD *) mycalloc (1, sizeof(particleMD));
//				AddItem1STD (&candidates, p1);			// register candidate (just one!)
				break;
			case LAYOUT_CYLINDER:
				p1 = (particleMD *) mycalloc (1, sizeof(particleMD));
				break;
			default:
				for (i=0; i<layoutList->n; i++) {		// peruse appropriate layout list
					p1 = layoutList->items[i].p1;		// pointer to particle
					if (p1->next)   continue;			// skip grafting sites
					if (p1->q != 0) continue;			// skip charged atom
					AddItem1STD (&candidates, p1);		// register atom
				}
		}

		// grow polymers
		for (i=0; i<polyM[set]; i++) {

			// grow polymer
			picked = 0;
			loop   = LOOP_MAX;
			while (!picked && loop--) {

				// no candidates left!
				if (candidates.n==0 && polyLayout[set]!=LAYOUT_ANCHOR && polyLayout[set]!=LAYOUT_FLUID && polyLayout[set]!=LAYOUT_PLATES && polyLayout[set]!=LAYOUT_CYLINDER && polyLayout[set]!=LAYOUT_RODX && polyLayout[set]!=LAYOUT_RODY && polyLayout[set]!=LAYOUT_U) error (EGRAFT);

				// choose candidate atom randomly
				c = (int) (RandomReal()*candidates.n);
				if (c==candidates.n) c--;
				if (c>0) p1 = candidates.items[c].p1;

				if (polyLayout[set] == LAYOUT_ANCHOR) {
					p1->rx = RandomReal()*sim->box[x_];
					p1->ry = RandomReal()*sim->box[y_];
					p1->rz = RandomReal()*sim->box[z_];
					p1->wx=p1->rx;
					p1->wy=p1->ry;
					p1->wz=p1->rz;
					p1->x0=p1->rx;
					p1->y0=p1->ry;
					p1->z0=p1->rz;
				}
				if (polyLayout[set] == LAYOUT_PLATES) {
					p1->rx = RandomReal()*sim->box[x_];
					p1->rz = RandomReal()*sim->box[z_];
					p1->wx = p1->rx;
					p1->wz = p1->rz;
					if (i%2==0) {
						p1->ry = 0+sim->sigma_lj;
						p1->wy = 0+sim->sigma_lj;
					}
					else {
						p1->ry = sim->box[y_]-sim->sigma_lj;
						p1->wy = sim->box[y_]-sim->sigma_lj;
					}

				}
				if (polyLayout[set]==LAYOUT_CYLINDER) {
					p1->rx = RandomReal()*sim->box[x_];
					theta = RandomReal()*2*M_PI;
					p1->ry = (sim->box[y_]-2)*sin(theta)/2.0+sim->box[y_]/2.0;
					p1->rz = (sim->box[z_]-2)*cos(theta)/2.0+sim->box[z_]/2.0;
					p1->wx = p1->rx;
					p1->wy = p1->ry;
					p1->wz = p1->rz;
					p1->mass = 1;
				}

				// assume we found a good site
				keep = 1;

				// skip site if it is too close to another polymer
				if (polySpread[set] > 0) {
					switch (polyLayout[set]) {

						case LAYOUT_SURFACE:
							// distance with other grafted polymers
							for (j=0; j<surface->n; j++) {
								p2 = surface->items[j].p1;
								if (!p2->next) continue;
								d2 = DistanceSquared (p1, p2, surfaceMetric, sim->box[y_]/2.0);
								if (d2 < d2Min) {
									keep = 0;
									break;
								}
							}
							break;

						case LAYOUT_FLUID:
							// distance with ALL other polymers
							// this may cause a crash...  I commented it out in another code
							for (j=0; j<polymer.n; j++) {
								p2 = polymer.items[j].p1;
								d2 = DistanceSquared (p1, p2, CARTESIAN, 0);
								if (d2 < d2Min) {
									keep = 0;
									break;
								}
							}
							break;

						// Tyler added - Identical to LAYOUT_FLUID (shouldn't be necessary)
						case LAYOUT_RODX:
							// distance with ALL other polymers
							// this may cause a crash...  I commented it out in another code
							for (j=0; j<polymer.n; j++) {
								p2 = polymer.items[j].p1;
								d2 = DistanceSquared (p1, p2, CARTESIAN, 0);
								if (d2 < d2Min) {
									keep = 0;
									break;
								}
							}
							break;

						// Zahra added - Identical to LAYOUT_RODX 
						case LAYOUT_RODY:
							// distance with ALL other polymers
							// this may cause a crash...  I commented it out in another code
							for (j=0; j<polymer.n; j++) {
								p2 = polymer.items[j].p1;
								d2 = DistanceSquared (p1, p2, CARTESIAN, 0);
								if (d2 < d2Min) {
									keep = 0;
									break;
								}
							}
							break;

						//Added by Karolina - identical to LAYOUT_FLUID
						case LAYOUT_U:
							// distance with ALL other polymers
							// this may cause a crash...  I commented it out in another code
							for (j=0; j<polymer.n; j++) {
								p2 = polymer.items[j].p1;
								d2 = DistanceSquared (p1, p2, CARTESIAN, 0);
								if (d2 < d2Min) {
									keep = 0;
									break;
								}
							}
							break;

						case LAYOUT_CYLINDER:
							// distance with ALL other polymers
							// this may cause a crash...  I commented it out in another code
							// there is something weird here, not sure why polymer.n-1 appears as the same atom
							// this function does not work, probably doesn't for LAYOUT_FLUID either...
							for (j=0; j<polymer.n; j++) {
								p2 = polymer.items[j].p1;
								if (!p2->next) continue;
								d2 = DistanceSquared (p1, p2, surfaceMetric, sim->box[y_]/2.0);
								if (d2 < d2Min) {
									keep = 0;
									break;
								}
							}
							break;

						case LAYOUT_PLATES:
							// distance with ALL other polymers
							// this may cause a crash...  I commented it out in another code
							// there is something weird here, not sure why polymer.n-1 appears as the same atom
							// this function does not work, probably doesn't for LAYOUT_FLUID either...
							for (j=0; j<polymer.n; j++) {
								p2 = polymer.items[j].p1;
								dx = p2->rx - p1->rx;
								dy = p2->ry - p1->ry;
								dz = p2->rz - p1->rz;
								ApplyPBC_dr (sim,&dx,&dy,&dz);
								dz = p2->rz - p1->rz; // I don't think this is necessary...
								d2 = dx*dx+dy*dy+dz*dz;
								if (d2 < d2Min) {
									keep = 0;
									break;
								}
							}
							break;

						case LAYOUT_ANCHOR:
							break;

						default:
							error (ELAYOUT);
					}
				}

				// grow it
				if (keep){
					if (polyLayout[set]==LAYOUT_FLUID ) {
						//this ia bad hack I think, it just avoids building from an existing atom
						//in my simulations I realized making multiple polymers cause cutting of bonds
						p1 = p3;
						p1->next = GrowLinearChain (sim, polyAtomType[set], layout, polyN[set],  NULL, &grown);
					}
					// Tyler added ---Stolen from above
					else if (polyLayout[set]==LAYOUT_RODX ) {
						p1 = p3;
						p1->next = GrowRodChain (sim, polyAtomType[set], layout, polyN[set],  NULL, &grown, x_);
					}
					// Zahra added ---Stolen from above
					else if (polyLayout[set]==LAYOUT_RODY ) {
						p1 = p3;
						p1->next = GrowRodChain (sim, polyAtomType[set], layout, polyN[set],  NULL, &grown, y_);
					}
					// added by Karolina - same as above
					else if (polyLayout[set]==LAYOUT_U ) {
						p1 = p3;
						p1->next = GrowUChain (sim, polyAtomType[set], layout, polyN[set],  NULL, &grown);
					}
					else p1->next = GrowLinearChain (sim, polyAtomType[set], layout, polyN[set], p1, &grown);
					if (grown) {
						// post-growth action depending on layout
						switch (polyLayout[set]) {
							case LAYOUT_SURFACE:
								// just register the polymer
								AddItemPoly (&polymer, p1, 1);
								break;
							case LAYOUT_FLUID:
								// detach from first fluid atom and register polymer
								// note this has been modded so we don't need to detach...
//								p2 = p1;
								p1 = p1->next;
								p1->prev = NULL;
//								p2->next = NULL;
								AddItemPoly (&polymer, p1, 0);
								break;
							// Tyler added --- Identical to LAYOUT_FLUID
							case LAYOUT_RODX:
								p1 = p1->next;
								p1->prev = NULL;
								AddItemPoly (&polymer, p1, 0);
								break;
							// Zahra added --- Identical to above
							case LAYOUT_RODY:
								p1 = p1->next;
								p1->prev = NULL;
								AddItemPoly (&polymer, p1, 0);
								break;
							//added by Karolina - identical to above
							case LAYOUT_U:
								p1 = p1->next;
								p1->prev = NULL;
								AddItemPoly (&polymer, p1, 0);
								break;
							case LAYOUT_PLATES:
								// detach from ghost atom, anchor and register polymer
								// for some reason using p2 instead of p3 below causes NAN in the FENE energies
								// perhaps because it wasn't calloc'ed?
								p3 = p1->next;
								p3->prev = NULL;
								p3->anchor = 1;
								p3->x0 = p1->wx;
								p3->y0 = p1->wy;
								p3->z0 = p1->wz;
								AddItemPoly (&polymer, p3, 1);
								break;
							case LAYOUT_CYLINDER:
								// detach from ghost atom, anchor and register polymer
								// for some reason using p2 instead of p3 below causes NAN in the FENE energies
								// perhaps because it wasn't calloc'ed?
								p3 = p1->next;
								p3->prev = NULL;
								p3->anchor = 1;
								p3->x0 = p1->wx;
								p3->y0 = p1->wy;
								p3->z0 = p1->wz;
								AddItemPoly (&polymer, p3, 1);
								break;
							case LAYOUT_ANCHOR:
								// detach from ghost atom, anchor and register polymer
								p2 = p1;
								p1 = p1->next;
								p1->prev = NULL;
								p1->anchor = 1;
								p1->x0 = p2->x0;
								p1->y0 = p2->y0;
								p1->z0 = p2->z0;
								AddItemPoly (&polymer, p1, 1);
								free(p2);
								break;
							default:
								error (ELAYOUT);
								break;
						}
						// note monomers are charged at the end of this function
						// we've picked a candidate
						picked = 1;
					}
				}

				// remove this candidate
				RemoveItem1STD (&candidates, c);
			}
			if (!picked) error (ELOOPMAX);
		}

		// report for this polymer set
		LOG ("  %-20s = %d\n", 	 			 "polyAtomType", 	polyAtomType[set]);
		LOG ("  %-20s = %d\n", 	 		 	 "polyLayout", 	 	polyLayout[set]);
		LOG ("  %-20s = %d\n", 	 			 "polyDensityKind", polyDensityKind[set]);
		LOG ("  %-20s = %f\n", 	 			 "polyDensity", 	polyDensity[set]);
		LOG ("  %-20s = %f\n", 	 			 "polySpread",	 	polySpread[set]);
		LOG ("  %-20s = %f\n", 	 			 "monoCharge",	 	monoCharge[set]);
		LOG ("  %-20s = %d polymers\n", 	 "polyM", 			polyM[set]);
		LOG ("  %-20s = %d monomers each\n", "polyN", 			polyN[set]);
		LOG ("  %-20s = %d monomers\n", 	 "Total", 			polyN[set]*polyM[set]);
	}

	// free temporary polymer list
	if (polymer.items)		free (polymer.items);
	if (candidates.items) 	free (candidates.items);
	ResetList1STD(&candidates);

	atom    = sim->atom.items;

	// report NB that only charge for 1 polymer set (set 0) is coded
	LOG ("Setting up charge on polymer\n");

	// rebuild candidate sites list (skipping busy atoms)
	ResetList1STD(&candidates);
	ResetList1STD(&charge);
	for (i=0; i<(sim->atom.n); i++) {		// peruse appropriate layout list
		p1 = atom + i;			// pointer to particle
		if (p1->q != 0) continue;	// skip charged atom
		if (p1->type != TYPE_MONOMER) continue;		// skip non monomers atom
		AddItem1STD (&candidates, p1);			// register atom
	}

	LOG ("%f = monocharge\n", monoCharge[0]);

	// put charges
	toAdd = fabs(monoCharge[0]*polyM[0]*polyN[0]);
//	toAdd = polyN[0]*polyM[0];
	loop  = -1;
	LOG ("%d = charges to add\n", toAdd);
	LOG ("%f = monoQ\n", monoCharge[0]);
	LOG ("%d = candidate monomers\n", candidates.n);
	while (toAdd && loop--) {
		if (candidates.n==0) error (ECHARGE);

		// choose candidate atom randomly
		c = (int) (RandomReal()*candidates.n);
		if (candidates.n==0) LOG ("out of candidates..\n");
		if (c==candidates.n) c--;
		p1 = candidates.items[c].p1;

		// assume success
		keep = 1;
		if (p1->type!=TYPE_MONOMER) keep = 0;
		// charge the atom
		if (keep) {
			AddItem1STD (&charge, p1);
			p1->q  = monoCharge[0]/fabs(monoCharge[0]);
//			p1->q  = monoCharge[0];
			toAdd--;
		}

		// remove this candidate
		RemoveItem1STD (&candidates, c);
		if (candidates.n==0) LOG ("something might still be wrong if charge !=/-1\n");
	}
	if (charge.items)		free (charge.items);
	if (candidates.items)	free (candidates.items);


}


/// Initializes electric charges. We attribute charges randomly to the
/// atoms. Note that if the fraction of charged atoms becomes large, this
/// method may become inefficient or fail altogether. All charges are given
/// in units of the elementary charge. We also try to avoid tight charge clusters
/// by trying to space out charges. Again, this is not efficient if there is
/// a large number of charges.
///
/// @param		sim	a pointer to a simulation structure
/// @return		void
/// @warning	This function may fail if the charge density is very high because
///				the atoms are picked at random. In such cases, the function
///				will fail and the program should exit with an ELOOPMAX.
/// @todo		Add an option to place charges according to a regular pattern.

//================================================================================
void InitCharges (simptr sim)
//================================================================================
{
	int			i, j, keep, set, c, toAdd, loop;
// 	int			layout;
	int			surfaceMetric;
	int			qBulkTot, qSurfaceTot;
	int			qNumber[Q_SETS], qCharge[Q_SETS], qLayout[Q_SETS], qDensityKind[Q_SETS];
	real		Qtot;
	real		qDensity[Q_SETS], qSpread[Q_SETS];
	real		A, V, d2, d2Min;
	particleMD	*p1, *p2;
	list1STD	*layoutList, *fluid, *surface, candidates;
	list1STD	charge;
	real		polyQCharge;

	// local sim variables
	A 	    = sim->surfaceArea;
	V 	    = sim->poreVolume;
	surface = &sim->surface;
	fluid   = &sim->fluid;

	// reset charge and candidate lists
	ResetList1STD(&charge);
	ResetList1STD(&candidates);

	// geometry variables
	surfaceMetric = CARTESIAN;
	switch (sim->geometry) {
		case GEOM_CAPILLARY:
			surfaceMetric = CYLINDRICAL | x_;
			break;
		default:
			break;
	}
// 	layout 	   = NONE;
	layoutList = NULL;
	d2Min 	   = 0;

	polyQCharge=0;
	for (i=0; i<sim->atom.n; i++) {
		polyQCharge += sim->atom.items[i].q;
	}

	polyQCharge = 0;
	LOG ("%f = polymer charge\n", polyQCharge);

	// charge set variables
	qBulkTot    = 0;
	qSurfaceTot = 0;
	for (set=0; set<Q_SETS; set++) {

		qLayout[set] = sim->qLayout[set];
		if (qLayout[set] == NONE) continue;

		qDensityKind[set] = sim->qDensityKind[set];
		qDensity[set] 	  = sim->qDensity[set];
		qSpread[set] 	  = sim->qSpread[set];
		qCharge[set]	  = sim->qCharge[set];
		qNumber[set]	  = sim->qNumber[set];

		switch (qDensityKind[set]) {
			case VOLUME:	qNumber[set] += qDensity[set] * V;
							break;
			case SURFACE:	qNumber[set] += qDensity[set] * A;
				if (qCharge[set] > 0 ) qNumber[set] -= polyQCharge;
				if (qNumber[set] < 0 ) {
					qCharge[set] = -qCharge[set];
					qNumber[set] = -qNumber[set];}
							break;
			default:		error (EDENSITY);
							break;
		}

		switch (qLayout[set]) {
			case LAYOUT_SURFACE:	qSurfaceTot+=qNumber[set];
									break;
			case LAYOUT_FLUID:		qBulkTot+=qNumber[set];
                       				break;
			// Tyler added --- Identical to LAYOUT_FLUID
			case LAYOUT_RODX:		qBulkTot+=qNumber[set];
                       				break;
			// Zahra added --- Identical to above
			case LAYOUT_RODY:		qBulkTot+=qNumber[set];
                       				break;
			// added by Karolina - identical to above
			case LAYOUT_U:		qBulkTot+=qNumber[set];
                       				break;
			default:			error (ELAYOUT);
							break;
		}
	}

	// build charge sets
	for (set=0; set<Q_SETS; set++) {

		// skip empty charge sets
		if (qNumber[set] <= 0 || qLayout[set]==NONE) continue;

		// report
		LOG ("Setting up charge set %d:\n", set);

		// setup some variables, according to the layout
		switch (qLayout[set]) {
			case LAYOUT_SURFACE:
				d2Min  	   = pow (qSpread[set]*sqrt(A/qSurfaceTot), 2.0);
// 				layout 	   = LAYOUT_SURFACE;
				layoutList = surface;
				break;
			case LAYOUT_FLUID:
				d2Min  	   = pow (qSpread[set]*pow(V/qBulkTot,1/3.0), 2.0);
// 				layout 	   = LAYOUT_FLUID;
				layoutList = fluid;
				break;
			// Tyler added --- Identical to LAYOUT_FLUID
			case LAYOUT_RODX:
				d2Min  	   = pow (qSpread[set]*pow(V/qBulkTot,1/3.0), 2.0);
// 				layout 	   = LAYOUT_FLUID;
				layoutList = fluid;
				break;
			// Zahra added --- Identical to above
			case LAYOUT_RODY:
				d2Min  	   = pow (qSpread[set]*pow(V/qBulkTot,1/3.0), 2.0);
// 				layout 	   = LAYOUT_FLUID;
				layoutList = fluid;
				break;
			// added by Karolina - identical to above
			case LAYOUT_U:
				d2Min  	   = pow (qSpread[set]*pow(V/qBulkTot,1/3.0), 2.0);
// 				layout 	   = LAYOUT_FLUID;
				layoutList = fluid;
				break;
			default:
				error (ELAYOUT)
				break;
		}

		// rebuild candidate sites list (skipping busy atoms)
		if (candidates.items) free (candidates.items);
		ResetList1STD(&candidates);
		for (i=0; i<layoutList->n; i++) {			// peruse appropriate layout list
			p1 = layoutList->items[i].p1;			// pointer to particle
			if (p1->next)   continue;				// skip grafting sites
			if (p1->q != 0) continue;				// skip charged atom
			AddItem1STD (&candidates, p1);			// register atom
		}

		// put charges
		toAdd = qNumber[set];
		loop  = LOOP_MAX;
		while (toAdd && loop--) {

			// choose candidate atom randomly
			c = (int) (RandomReal()*candidates.n);
			if (c==candidates.n) c--;
			p1 = candidates.items[c].p1;

			// assume success
			keep = 1;

			// skip when too close to another charge (according to geometry)
			if (qSpread[set] > 0) {
				switch (qLayout[set]) {

					case LAYOUT_SURFACE:
						// distance to other surface charges
						for (j=0; j<surface->n; j++) {
							p2 = surface->items[j].p1;
							if (p2->q == 0) continue;
							d2 = DistanceSquared (p1, p2, surfaceMetric, sim->box[y_]/2.0);
							if (d2 < d2Min) {
								keep = 0;
								break;
							}
						}
						break;

					case LAYOUT_FLUID:
						// distance with ALL other charges
						for (j=0; j<charge.n; j++) {
							p2 = charge.items[j].p1;
							d2 = DistanceSquared (p1, p2, CARTESIAN,0);
							if (d2 < d2Min) {
								keep = 0;
								break;
							}
						}
						break;

					// Tyler added --- Identical to LAYOUT_FLUID
					case LAYOUT_RODX:
						// distance with ALL other charges
						for (j=0; j<charge.n; j++) {
							p2 = charge.items[j].p1;
							d2 = DistanceSquared (p1, p2, CARTESIAN,0);
							if (d2 < d2Min) {
								keep = 0;
								break;
							}
						}
						break;

					// Zahra added --- Identical to above
					case LAYOUT_RODY:
						// distance with ALL other charges
						for (j=0; j<charge.n; j++) {
							p2 = charge.items[j].p1;
							d2 = DistanceSquared (p1, p2, CARTESIAN,0);
							if (d2 < d2Min) {
								keep = 0;
								break;
							}
						}
						break;

					// added by Karolina - identical to above
					case LAYOUT_U:
						// distance with ALL other charges
						for (j=0; j<charge.n; j++) {
							p2 = charge.items[j].p1;
							d2 = DistanceSquared (p1, p2, CARTESIAN,0);
							if (d2 < d2Min) {
								keep = 0;
								break;
							}
						}
						break;
					default:
						error (ELAYOUT);
				}
			}

			// charge the atom
			if (keep) {
				AddItem1STD (&charge, p1);
				p1->q  = qCharge[set];
				toAdd--;
			}

			// remove this candidate
			RemoveItem1STD (&candidates, c);
			if (candidates.n==0) error (ECHARGE);
		}

		// report for this charge set
		LOG ("  %-20s = %d\n", 	 "qLayout", 	 	qLayout[set]);
		LOG ("  %-20s = %d\n", 	 "qDensityKind", 	qDensityKind[set]);
		LOG ("  %-20s = %f\n", 	 "qDensity", 		qDensity[set]);
		LOG ("  %-20s = %f\n", 	 "qSpread",	 		qSpread[set]);
		LOG ("  %-20s = %d\n", 	 "qCharge",	 		qCharge[set]);
		LOG ("  %-20s = %d\n",	 "qNumber", 		qNumber[set]);
	}

	// free temporary charge list
	if (charge.items)		free (charge.items);
	if (candidates.items)	free (candidates.items);

	// report overall charge
	Qtot=0;
	for (i=0; i<sim->atom.n; i++) {
		Qtot += sim->atom.items[i].q;
	}
	LOG ("Charge initialization complete\n");
	LOG ("  %-20s = %f\n",	 "Net Charge", Qtot);
	if (Qtot != 0) LOG ("*** WARNING *** The system has a non-zero net charge!\n");
	if (Qtot != 0) LOG ("*** WARNING *** The system has a non-zero net charge!\n");
	if (fabs(Qtot)>1) {
		LOG ("You have some serious net charge issues, check your parameters.\n");
//		exit(0);
	}
	LOG ("Successfully built charges\n");

}


/// Initializes the velocities of all the particles. Each particle is given
/// the mean velocity (set by the temperature) in a random direction. There
/// is 3*n-3 degrees of freedom (because the total momentum is conserved) so
/// from the equipartition theorem we can equate 3*(n-1)*kT = n*v^2
/// [Rapaport]. Actually, if the temperature is kept constant then there is
/// an additional constraint, so we really have (3*n-4) degrees of freedom
/// [Heermann]; the discrepancy is irrelevant when nAtom is large of course.
///
/// @param		sim a pointer to a simulation structure
/// @return		void

//================================================================================
void InitVelocities (simptr sim)
//================================================================================
{
	int			i, nAtom;
	particleMD	*atom;
	real		kT, kinEMean, vMag=0;
	real		vector[DIM_MD], Vcm[DIM_MD], Mtot;

	// local sim variables
	atom  = sim->atom.items;
	nAtom = sim->atom.n;
	kT	  = sim->kT[sim->phase];

	// report
	LOG ("Initializing particle velocities\n");

	// compute mean velocity from temperature (3*nAtom-3 degrees of freedom)
	kinEMean = 0.5 * DIM_MD * (1-1.0/nAtom) * kT;

	// set the center-of-mass velocity to zero
	Vcm[x_] = 0;
	Vcm[y_] = 0;
	Vcm[z_] = 0;
	Mtot	= 0;

	// set randomly oriented velocities
	for (i=0; i<nAtom; i++) {
		RandomVector3D (vector);
		vMag = sqrt (2.0 * kinEMean / atom[i].mass);
		atom[i].vx = vMag * vector[x_];
		atom[i].vy = vMag * vector[y_];
		atom[i].vz = vMag * vector[z_];
		Vcm[x_] += atom[i].vx * atom[i].mass;
		Vcm[y_] += atom[i].vy * atom[i].mass;
		Vcm[z_] += atom[i].vz * atom[i].mass;
		Mtot	+= atom[i].mass;
	}

	// adjust velocities so that CM is at rest
	Vcm[x_] /= Mtot;
	Vcm[y_] /= Mtot;
	Vcm[z_] /= Mtot;
	for (i=0; i<nAtom; i++) {
		atom[i].vx -= Vcm[x_];
		atom[i].vy -= Vcm[y_];
		atom[i].vz -= Vcm[z_];
	}

	// report
	LOG ("  %-20s = %g\n", "kinEMean", kinEMean);
	LOG ("  %-20s = (%g, %g, %g)\n", "Vcm shift", -Vcm[x_], -Vcm[y_], -Vcm[z_]);
}


/// Builds a list of pointers to particles that constitute different parts
/// of the system. This list is mostly useful in setting up particles in a specific
/// region of the system.
///
///	@param		sim a pointer to a simulation structure
/// @return		void
/// @todo		Add the option of defining many independent regions sets in the
///				simulation, in a more robust way.

//================================================================================
void SetupLayoutLists (simptr sim)
//================================================================================
{
	int 		i, nAtom;
	particleMD	*atom, *p;
	list1STD	*wall, *fluid, *surface;
	real		r2, rMax2;

	// local sim variables
	atom	 = sim->atom.items;
	nAtom	 = sim->atom.n;
	wall     = &sim->wall;
	fluid    = &sim->fluid;
	surface	 = &sim->surface;

	// reset surface list
	ResetList1STD (wall);
	ResetList1STD (fluid);
	ResetList1STD (surface);

	// build surface list, according to geometry
	switch (sim->geometry) {

		case GEOM_BULK:
			// loop over all atoms
			for (i=0; i<nAtom; i++) {
				p = atom+i;
				if (p->type == TYPE_WALL)				AddItem1STD (wall, p);
				if (p->type == TYPE_FLUID)				AddItem1STD (fluid, p);
			}
			break;

		case GEOM_PLATES:
			// loop over all atoms
			for (i=0; i<nAtom; i++) {
				p = atom+i;
				if (p->type == TYPE_WALL)				AddItem1STD (wall, p);
				if (p->type == TYPE_FLUID)				AddItem1STD (fluid, p);
			}
			break;
		case GEOM_CYLINDER:
			// loop over all atoms
			for (i=0; i<nAtom; i++) {
				p = atom+i;
				if (p->type == TYPE_WALL)				AddItem1STD (wall, p);
				if (p->type == TYPE_FLUID)				AddItem1STD (fluid, p);
			}
			break;

		case GEOM_CAPILLARY:
			// maximum radius for inside capillary surface
			rMax2 = pow (sim->caprInTrue+0.75,2.0);
			// loop over all atoms
			for (i=0; i<nAtom; i++) {
				p = atom+i;
				r2 = p->ry*p->ry + p->rz*p->rz;
				if (p->type == TYPE_WALL)				AddItem1STD (wall, p);
				if (p->type == TYPE_FLUID)				AddItem1STD (fluid, p);
				if (p->type == TYPE_WALL && r2 < rMax2)	AddItem1STD (surface, p);
			}
			break;

		default:
			error(EGEOMETRY);
			break;
	}
}


//================================================================================
void SetupPolymerList (simptr sim)
//================================================================================
{
	int 		i, nAtom, grafted;
	particleMD	*atom, *p1;

	// local sim variables
	atom  = sim->atom.items;
	nAtom = sim->atom.n;

	// report
	LOG ("Building list of polymers\n");

	// reset simulation charge list
	ResetListPoly (&sim->polymer);

	// build list of polymers
	for (i=0; i<nAtom; i++) {
		p1 = atom+i;
		grafted = 0;
		if (p1->type == TYPE_MONOMER && p1->prev == NULL) {
			if (p1->anchor == 1) grafted = 1;
			AddItemPoly (&sim->polymer, p1, grafted);
		}
	}

	// report
	LOG ("  %-20s = %d\n", "nPolymer", sim->polymer.n);
}


//================================================================================
void SetupChargeList (simptr sim)
//================================================================================
{
	// Builds a list of all charges in the simulation. The list is an array of
	// pointers to all particles for which q != 0. This list is used in
	// computing electrostatic interactions in the main force loop.

	int 		i, nAtom;
	particleMD	*atom, *p1;

	// local sim variables
	atom  = sim->atom.items;
	nAtom = sim->atom.n;

	// report
	LOG ("Building list of charges\n");

	// reset simulation charge list
	ResetList1STD (&sim->charge);

	// build list of charges
	for (i=0; i<nAtom; i++) {
		p1 = atom+i;
		if (p1->q != 0) AddItem1STD (&sim->charge, p1);
	}

	// report
	LOG ("  %-20s = %d\n", "nCharge", sim->charge.n);
}


//================================================================================
void SetupNeighborList (simptr sim)
//================================================================================
{
	// Initialize the neighbor management system of the simulation. Calls other
	// functions to allocate the cell list and the pbc list, and to refresh the
	// neighbor list (i.e. build it for the first time).

	int		i, nCell, *nCellAxis;
// 	int		nAtom;
	real	rCut, rNebrShell, rNebr, *box, *cellInvWidth;

	// local sim variables
// 	nAtom	 	 = sim->atom.n;
	rCut 	 	 = sim->rCut;
	rNebrShell	 = sim->rNebrShell;
	nCellAxis 	 = sim->nCellAxis;
	cellInvWidth = sim->cellInvWidth;
	box		 	 = sim->box;

	// report
	LOG ("Building list of short-range neighbors\n");

	// reset neighbor lists
	ResetList2STD (&sim->nebrSTD);
	ResetList2PBC (&sim->nebrPBC);

	// calculate the neighborhood radius
	rNebr = (rCut + rNebrShell);
	sim->rNebr  = rNebr;
	sim->rNebr2 = rNebr*rNebr;

	// calculate the number and size of cells
	nCell = 1;
	for (i=0; i<DIM_MD; i++) {
		nCellAxis[i] = box[i] / rNebr;
		nCell *= nCellAxis[i];
		cellInvWidth[i] = nCellAxis[i]/box[i];
	}
	sim->nCell = nCell;

	// allocate the cell and PBC offset lists
	AllocateCellList (sim);
	AllocatePBCList  (sim);

	// build the neighbor list
	RefreshNebrList (sim);

	// report
	LOG ("  %-20s = %g\n", "rNebr", rNebr);
	LOG ("  %-20s = %d\n", "nCell", nCell);
	LOG ("  %-20s = %d\n", "nNebrSTD", sim->nebrSTD.n);
	LOG ("  %-20s = %d\n", "nNebrPBC", sim->nebrPBC.n);
	LOG ("  %-20s = %d\n", "nNebrTOT", sim->nebrSTD.n + sim->nebrPBC.n);
}


//================================================================================
void SetupAnchorList (simptr sim)
//================================================================================
{
	// Initializes a list of the atoms that are attached to a fixed crystal
	// lattice throughout the simulation (usually these are the wall atoms).

	int			n, nAtom;
	particleMD	*atom;

	// local sim variables
	atom  = sim->atom.items;
	nAtom = sim->atom.n;

	// report
	LOG ("Building list of achored atoms\n");

	// reset anchor list
	ResetList1STD (&sim->anchor);

	// build list of all anchored atom pointers
	for (n=0; n<nAtom; n++) {
		if (atom[n].anchor) AddItem1STD (&sim->anchor, atom+n);
	}
	// report
	LOG ("  %-20s = %d\n", "nAnchor", sim->anchor.n);
}


//================================================================================
void SetupFeneList (simptr sim)
//================================================================================
{
	// Initializes a list of pairs that are bounded by a FENE interaction.

	int			n, nn, i, nAtom, nPolymer, nMonomer;
	particleMD	*atom, *p;

	// local sim variables
	atom  = sim->atom.items;
	nAtom = sim->atom.n;
	// need to consider multiple polymer chains
	nPolymer = sim->polyM[0];
	nMonomer = sim->polyN[0];

	// report
	LOG ("Building list of FENE pairs\n");

	// reset fene pairs list
	ResetList2STD (&sim->fene);

	// build list of all anchored atom pointers
	for(i=0 ; i<nPolymer ;i++){
		nn = i*nMonomer ;
		for (n=nn; n<nn+nMonomer-1; n++) {
			p = atom+n;
			AddItem2STD (&sim->fene, p, p->next);
		}
	}

	// report
	LOG ("  %-20s = %d\n", "nFene", sim->fene.n);
}


//================================================================================
void SetupBendList (simptr sim)
//================================================================================
{
	// Initializes a list of pairs that are bounded by a bend interaction.

	int			n, nn, i, nAtom, nPolymer, nMonomer;
	particleMD	*atom, *p;

	// local sim variables
	atom  = sim->atom.items;
	nAtom = sim->atom.n;
	// need to consider multiple polymer chains
	nPolymer = sim->polyM[0];
	nMonomer = sim->polyN[0];

	// report
	LOG ("Building list of bend triplets\n");

	// reset bend triplets list
	ResetList3STD (&sim->bend);

	// build list of all anchored atom pointers
	if (nAtom>=3){
		for(i=0 ; i<nPolymer ;i++){
			nn = 1 + i*nMonomer ;
			for (n=nn; n<nn+nMonomer-2; n++) {
				p = atom+n;
				AddItem3STD (&sim->bend, p->prev, p, p->next);
			}
		}
	}

	// report
	LOG ("  %-20s = %d\n", "nBend", sim->bend.n);
}


//================================================================================
void SetupGroups (simptr sim)
//================================================================================
{
	// Assign group numbers to groups of atoms. Groups are defined in the header
	// file "mdtypes.h". Group numbers consist of a bitfield so that atoms
	// can concisely be part of many groups simultaneously.

	int			n, nAtom, nPolymer, group, count[GROUP_MAX];
	char		groupStr[STRLEN];
	real		y, z, r, rmid;
	particleMD	*atom, *p;
	itemPoly	*polymer;

	// local sim variables
	atom  	 = sim->atom.items;
	nAtom 	 = sim->atom.n;
	polymer  = sim->polymer.items;
	nPolymer = sim->polymer.n;
	rmid  	 = sim->caprIn + 0.5*(sim->caprOut-sim->caprIn);

	// report
	LOG ("Setting up atom groups\n");

	// reset group counters
	for (n=0; n<GROUP_MAX; n++) count[n] = 0;

	// set atom groups
	for (n=0; n<nAtom; n++) {

		// reset group
		atom[n].group = GROUP_NONE;

		// fluid group
		if (atom[n].type  == TYPE_FLUID) {
			atom[n].group |= GROUP_FLUID;
			count[GroupToIndex(GROUP_FLUID)]++;
		}

		// wall groups
		if (atom[n].type  == TYPE_WALL) {
			atom[n].group |= GROUP_WALL;
			count[GroupToIndex(GROUP_WALL)]++;
			y = atom[n].ry;
			z = atom[n].rz;
			r = sqrt (y*y + z*z);
			if (r < rmid) {
				atom[n].group |= GROUP_WALL_IN;
				count[GroupToIndex(GROUP_WALL_IN)]++;
			}
			else {
				atom[n].group |= GROUP_WALL_OUT;
				count[GroupToIndex(GROUP_WALL_OUT)]++;
			}
		}

		// monomer group
		if (atom[n].type  == TYPE_MONOMER) {
			atom[n].group |= GROUP_MONOMER;
			count[GroupToIndex(GROUP_MONOMER)]++;
		}

		// positive ions
		if (atom[n].q > 0) {
			atom[n].group |= (GROUP_ION | GROUP_ION_POS);
			count[GroupToIndex(GROUP_ION_POS)]++;
			count[GroupToIndex(GROUP_ION)]++;
		}

		// negative ions
		if (atom[n].q < 0) {
			atom[n].group |= (GROUP_ION | GROUP_ION_NEG);
			count[GroupToIndex(GROUP_ION_NEG)]++;
			count[GroupToIndex(GROUP_ION)]++;
		}
	}

	// set object types (for visual appearance in viewer)
	for (n=0; n<nAtom; n++) {
		if		(atom[n].group & GROUP_ION_POS)	atom[n].object = OBJ_ION_POS;
		else if (atom[n].group & GROUP_ION_NEG)	atom[n].object = OBJ_ION_NEG;
		else if (atom[n].group & GROUP_WALL)	atom[n].object = OBJ_WALL;
		else if (atom[n].group & GROUP_FLUID)	atom[n].object = OBJ_FLUID;
		else if (atom[n].group & GROUP_MONOMER)	atom[n].object = OBJ_MONOMER;
	}

	// set object types for first monomers (all)
	for (n=0; n<nPolymer; n++) {
		if (polymer[n].p1) polymer[n].p1->object = OBJ_MONOMER_1;
	}

	// change object types for grafted polymers
	for (n=0; n<nPolymer; n++) {
		if (polymer[n].grafted) {
			p = polymer[n].p1;
			p->object = OBJ_MONOMER_GRAFT_1;
			p->group |= GROUP_GRAFT;
			count[GroupToIndex(GROUP_GRAFT)]++;
			p = p->next;
			while (p) {
				p->object = OBJ_MONOMER_GRAFT;
				p = p->next;
			}
		}
	}

	// report
	for (n=1; n<GROUP_MAX; n++) {
		if (count[n] > 0) {
			group = GroupFromIndex(n);
			GroupToStr (group, groupStr);
			LOG ("  %-14s (%#010X) = %d atoms\n", groupStr, group, count[n]);
		}
	}
}


//================================================================================
void SetupMeasurements (simptr sim)
//================================================================================
{
	// Initializes variables related to the data collection process.

	// reset properties accumulators
	AccumProperties (sim, ACTION_RESET);

	// setup measurement related lists
	SetupHistograms (sim);
// 	SetupScenes		(sim);
	SetupVMD		(sim);
}


//================================================================================
void SetupStepCounters (simptr sim)
//================================================================================
{
	// Allocates memory and initializes the array of pointers to all step
	// counters of the simulation. This function needs to be modified when a new
	// step counter is introduced in the simulation. Simply increment add
	// the corresponding NewCounter call.

	// reset counters
	sim->stepCounter.items = NULL;
	sim->stepCounter.n     = 0;
	sim->stepCounter.max   = 0;

	// global step counter (do not remove);
	AddItemStep (&sim->stepCounter, sim->step, NULL);

	// other counters
	AddItemStep (&sim->stepCounter, sim->stepAvg, 		   &AverageProperties);
	AddItemStep (&sim->stepCounter, sim->stepCheckpoint,   &CheckpointWrite);
	AddItemStep (&sim->stepCounter, sim->stepThermRescale, &ThermostatRescale);
}


/// Slowly relax the system to make sure atoms are not too close to each
/// other. The function also prints the initial scene as well as the relaxed state.
/// This function only contains the iteration loop. The actual relaxation step
/// is carried out by the RelaxStep function.
///
/// @param		sim a pointer to a simulation structure
/// @return		void

//================================================================================
void RelaxParticles (simptr sim)
//================================================================================
{
	int 	 stepRelax, mdstep;
	sceneptr s;

	int	 	i, nAtom;
	particleMD	*atom, *p;
	real		dt, v2, v2max;
	real maxFene2, bx, by, bz, b2;

	// local sim variables
	atom	= sim->atom.items;
	nAtom	= sim->atom.n;
	dt		= sim->dt;

	// local sim variables
	stepRelax = sim->stepRelax;
	maxFene2 = 0.9*sim->r0Fene * 0.9*sim->r0Fene;

	// report
	LOG ("Relaxation of particle coordinates\n");
	LOG ("  %-20s = %d\n", "stepRelax", sim->stepRelax);
	LOG ("  %-20s = %f\n", "drRelax", 	sim->drRelax);

	// get relaxation scene
// 	s = SceneGet (sim->scenes, "scene-relax");
	s = VMDGet (sim->scenes, "vmd-relax");

	// print the initial state
// 	SceneFunction ((void *) sim, s, ACTION_PRINT);
	VMDFunction ((void *) sim, s, ACTION_PRINT);

// 	for (i=0; i<nAtom; i++) {
// 		p = atom+i;
// 		p->rx = i;
// 		p->ry = 5.;
// 		p->rz = 5.;
// 		p->wx = i;
// 		p->wy = 5.;
// 		p->wz = 5.;
// 		p->vx = 0.01;
// 		p->vy = 0.01;
// 		p->vz = 0.01;
// 	}

	// relaxation loop
/*
	while (stepRelax--) {

		// take one relaxation step
		RelaxStep (sim);

		// update neighbor list if needed
		CheckNebrList (sim);

		// feedback
		LOG ("%d\n", stepRelax);
	}
*/
	stepRelax = sim->stepRelax;

	// report
	LOG ("Relaxation of particle coordinates using a force cap...\n");
	LOG ("  %-20s = %d\n", "ForceCapRelaxSteps", sim->stepRelax);

	while (stepRelax--) {

		for (mdstep=0;mdstep<10;mdstep++) {
			v2max = 0;

			ComputeCapForces (sim);

			for (i=0; i<nAtom; i++) {
				p = atom+i;
				if (p->type==TYPE_WALL) continue;

				// enforce maximum FENE bond lengths
				if (p->next) {
					bx = p->wx+(p->vx+p->ax*dt)*dt - p->next->wx-(p->next->vx+p->next->ax*dt)*dt;
					by = p->wy+(p->vy+p->ay*dt)*dt - p->next->wy-(p->next->vy+p->next->ay*dt)*dt;
					bz = p->wz+(p->vz+p->az*dt)*dt - p->next->wz-(p->next->vz+p->next->az*dt)*dt;
					b2 = bx*bx + by*by + bz*bz;
					if (b2 > maxFene2) continue;
				}
				if (p->prev) {
					bx = p->wx+(p->vx+p->ax*dt)*dt - p->prev->wx-(p->prev->vx+p->prev->ax*dt)*dt;
					by = p->wy+(p->vy+p->ay*dt)*dt - p->prev->wy-(p->prev->vy+p->prev->ay*dt)*dt;
					bz = p->wz+(p->vz+p->az*dt)*dt - p->prev->wz-(p->prev->vz+p->prev->az*dt)*dt;
					b2 = bx*bx + by*by + bz*bz;
					if (b2 > maxFene2) continue;
				}
				if (p->next) {
					bx = p->wx - p->next->wx;
					by = p->wy - p->next->wy;
					bz = p->wz - p->next->wz;
					b2 = bx*bx + by*by + bz*bz;
					if (b2 > maxFene2) exit(0);
				}
				if (p->prev) {
					bx = p->wx - p->prev->wx;
					by = p->wy - p->prev->wy;
					bz = p->wz - p->prev->wz;
					b2 = bx*bx + by*by + bz*bz;
					if (b2 > maxFene2) exit(0);
				}

				// increment velocities
				p->vx = p->ax*dt;
				p->vy = p->ay*dt;
				p->vz = p->az*dt;
				// increment positions
				p->rx += p->vx*dt;
				p->ry += p->vy*dt;
				p->rz += p->vz*dt;
				// increment world (true) positions
				p->wx += p->vx*dt;
				p->wy += p->vy*dt;
				p->wz += p->vz*dt;
				v2 = (p->vx*p->vx) + (p->vy*p->vy) + (p->vz*p->vz);
				if (v2 > v2max) v2max = v2;
//				if (v2 > 0.01*sim->kT[1]) mdstep = 0;
			}


			// update neighbor list if needed
			sim->drTotMax += sqrt(v2max)*dt;
			if (v2max>10) {
				LOG ("what is wrong??? %E\n", v2max);
			}
// 			CheckNebrList (sim);
 			PeriodicBoundaries (sim);
		}


		// feedback
		LOG ("%d\n", stepRelax);
		sim->stepRelax = stepRelax;

	}
	sim->stepRelax = 1;
	for (mdstep=0;mdstep<0;mdstep++) {
		v2max = 0;
		ComputeCapForces (sim);
// 		ComputeForces (sim);
		LOG ("setting reasonable velocities for thermostat %d\n", mdstep);
		for (i=0; i<nAtom; i++) {
			p = atom+i;
//			if (p->type==TYPE_WALL) continue;
			p->vx = p->ax*dt;
			p->vy = p->ay*dt;
			p->vz = p->az*dt;
			if (p->next) {
				bx = p->wx+(p->vx+p->ax*dt)*dt - p->next->wx-(p->next->vx+p->next->ax*dt)*dt;
				by = p->wy+(p->vy+p->ay*dt)*dt - p->next->wy-(p->next->vy+p->next->ay*dt)*dt;
				bz = p->wz+(p->vz+p->az*dt)*dt - p->next->wz-(p->next->vz+p->next->az*dt)*dt;
				b2 = bx*bx + by*by + bz*bz;
				if (b2 > maxFene2) continue;
			}
			if (p->prev) {
				bx = p->wx+(p->vx+p->ax*dt)*dt - p->prev->wx-(p->prev->vx+p->prev->ax*dt)*dt;
				by = p->wy+(p->vy+p->ay*dt)*dt - p->prev->wy-(p->prev->vy+p->prev->ay*dt)*dt;
				bz = p->wz+(p->vz+p->az*dt)*dt - p->prev->wz-(p->prev->vz+p->prev->az*dt)*dt;
				b2 = bx*bx + by*by + bz*bz;
				if (b2 > maxFene2) continue;
			}
			if (p->next) {
				bx = p->wx - p->next->wx;
				by = p->wy - p->next->wy;
				bz = p->wz - p->next->wz;
				b2 = bx*bx + by*by + bz*bz;
				if (b2 > maxFene2) exit(0);
			}
			if (p->prev) {
				bx = p->wx - p->prev->wx;
				by = p->wy - p->prev->wy;
				bz = p->wz - p->prev->wz;
				b2 = bx*bx + by*by + bz*bz;
				if (b2 > maxFene2) exit(0);
			}
			// increment positions
			p->rx += p->vx*dt;
			p->ry += p->vy*dt;
			p->rz += p->vz*dt;
			// increment world (true) positions
			p->wx += p->vx*dt;
			p->wy += p->vy*dt;
			p->wz += p->vz*dt;
			v2 = (p->vx*p->vx) + (p->vy*p->vy) + (p->vz*p->vz);
			if (v2 > v2max)
				v2max = v2;
		}
		// update neighbor list if needed
		sim->drTotMax += sqrt(v2max)*dt;
		if (v2max>10) {
			LOG ("what is really wrong??? %E\n", v2max);
		}
// 		CheckNebrList (sim);
		PeriodicBoundaries (sim);
	}

	PeriodicBoundaries (sim);

 	for (i=0; i<nAtom; i++) {
 		p = atom+i;
// 		p->rx = i;
// 		p->ry = 5.;
// 		p->rz = 5.;
// 		p->wx = i;
// 		p->wy = 5.;
// 		p->wz = 5.;
 		p->vx = 0.0;
 		p->vy = 0.0;
 		p->vz = 0.0;
 	}

// 	RefreshNebrList (sim);

}


/// Move each atom a small distance in the direction of the force, in order
/// to relax the system in a state appropriate for MD integration. The usual
/// ComputeForces function is used, but instead of the standard MD integration
/// loop, we simply move each atom by a small predefined distance (controlled
/// by the input parameter drRelax) in the direction of the net force. In effect,
/// this is a cheap steepest descent method!
///
/// @param		sim a pointer to a simulation structure
/// @return		void
/// @warning	If atoms are too close initially, the MD integration may "explode".
///				This is because is the acceleration is too great, the cpu may
///				return 1/Inf ~ 0 for the length of the vector, without warning.
///				Hence, such atom pairs are not separated and the standard MD loop
///				then crashes within the first few steps of the simulation.

//================================================================================
void RelaxStep (simptr sim)
//================================================================================
{
	int			i, nAtom;
	particleMD 	*atom, *p;
	real		a,  ax, ay, az;
	real		dx, dy, dz;
	real		b2, bx, by, bz;
	real		drRelax, maxFene2;
	real		scale;

	// local sim variables
	atom	 = sim->atom.items;
	nAtom	 = sim->atom.n;
	drRelax  = sim->drRelax;
	maxFene2 = 0.9*sim->r0Fene * 0.9*sim->r0Fene;

	// calculate the net force on all particles
	ComputeForces (sim);

	// move each atom in direction of net force
	for (i=0; i<nAtom; i++) {

		// get accelerations
		p  = atom+i;
		ax = p->ax;
		ay = p->ay;
		az = p->az;

		// scale down acceleration (to avoid floating point errors!)
		scale = fabs(ax) + fabs(ay) + fabs(az);
		if (scale > 10) {
			ax /= scale;
			ay /= scale;
			az /= scale;
		}
		a  = ax*ax + ay*ay + az*az;

		// calculate displacement
		if (a > 0) {

			// acceleration vector norm
			a = sqrt(a);

			// unit vector in direction of net force
			ax /= a;
			ay /= a;
			az /= a;

			// small relaxation displacement
			dx = ax * drRelax;
			dy = ay * drRelax;
			dz = az * drRelax;

			// enforce maximum FENE bond lengths
			if (p->next) {
				bx = p->wx+dx - p->next->wx;
				by = p->wy+dy - p->next->wy;
				bz = p->wz+dz - p->next->wz;
				b2 = bx*bx + by*by + bz*bz;
				if (b2 > maxFene2) continue;
			}
			if (p->prev) {
				bx = p->wx+dx - p->prev->wx;
				by = p->wy+dy - p->prev->wy;
				bz = p->wz+dz - p->prev->wz;
				b2 = bx*bx + by*by + bz*bz;
				if (b2 > maxFene2) continue;
			}

			// move atom
			p->rx += dx;
			p->ry += dy;
			p->rz += dz;
			p->wx += dx;
			p->wy += dy;
			p->wz += dz;
		}
	}

	// update maximum displacement
	sim->drTotMax += sim->drRelax;
}


/// Account for the possibility that particles may overlap within the
/// distance drMin (should be ~ 0.01), and if so displace one of them
/// slightly to avoid divergences in the initial integration steps.
/// Particularly important when using a float incarnation of the code.
///
/// @param		sim a pointer to a simulation structure
/// @return		void

//================================================================================
void JiggleParticles (simptr sim)
//================================================================================
{
    int         i, n, overlap;
    real     	dx, dy, dz, dr2, drMin, drMin2;
    real 	    *pbc;
    particleMD    *atom, *p1, *p2;
    item2STD    *nebrSTD;
    item2PBC    *nebrPBC;

    // local sim variables
    atom     = sim->atom.items;
    nebrSTD  = sim->nebrSTD.items;
    nebrPBC  = sim->nebrPBC.items;
    drMin    = sim->overlapMin;
    drMin2   = drMin*drMin;

    // report
    LOG ("Jiggling particles to remove overlaps\n");

    // check overlap
    do {

		// assume there is no overlap
        overlap=0;

		// check STD neigbors
        n = sim->nebrSTD.n;
        for (i=0; i<n; i++) {

            // extract particle pointers
            p1 = nebrSTD[i].p1;         		// pointer to the first particle
            p2 = nebrSTD[i].p2;         		// pointer to the second particle

            // calculate dr
            dx  = (p2->rx - p1->rx);
            dy  = (p2->ry - p1->ry);
            dz  = (p2->rz - p1->rz);
            dr2 = (dx*dx + dy*dy + dz*dz);

            // check overlap and jiggle if needed
            if (dr2 < drMin2) {
                overlap = 1;
                p2->rx += 2*drMin;
				p2->wx += 2*drMin;
				LOG ("  jiggled particle %d: rx += %f\n", (int)(p2-atom), 2*drMin);
            }
        }

		// check PBC neighbors
        n = sim->nebrPBC.n;
        for (i=0; i<n; i++) {

            // extract particle pointers
		p1  = nebrPBC[i].p1;                // pointer to the first particle
            p2  = nebrPBC[i].p2;                // pointer to the second particle
            pbc = nebrPBC[i].pbc;               // pointer to the PBC offset

            // calculate dr
            dx = (p2->rx - p1->rx) + pbc[x_];
            dy = (p2->ry - p1->ry) + pbc[y_];
            dz = (p2->rz - p1->rz) + pbc[z_];
            dr2 = (dx*dx + dy*dy + dz*dz);

            // check overlap and jiggle if needed
            if (dr2 < drMin2) {
                overlap = 1;
                p2->rx += 2*drMin;
				p2->wx += 2*drMin;
				LOG ("  jiggled particle %d: rx += %f\n", (int)(p2-atom), 2*drMin);
            }
        }

		// update maximum displacement and check neighbor lists
		sim->drTotMax += 2*drMin;
// 		CheckNebrList(sim);
		PeriodicBoundaries (sim);

    } while (overlap);
}


/// Housekeeping function at the end of the simulation. This function is
/// called just after the main md simulation loop, so put wrap-up code in here.
///
///	@param		sim a pointer to a simulation structure
///	@return		void

//================================================================================
void FinishSimulation (simptr sim)
//================================================================================
{
	// report
	LOG ("--------------------------------------------------------------\n");
	LOG ("END of main simulation loop\n");
	LOG ("--------------------------------------------------------------\n");

	// final calculations
	FinalCalculations (sim);

	// close data files
	CloseDataFiles (sim);

	// release dynamically allocated memory
	FreeMemory (sim);

	// remove pid file
	unlink ("pid");

	// terminate MPI execution
	#ifdef MPI
	MPI_Finalize();
	#endif
}


//================================================================================
particleMD *GrowLinearChain (simptr sim, int type, int layout, int n, particleMD *p0, int *status)
//================================================================================
{
	// RECURSIVE function to grow a n-monomer subchain from the particle pointed
	// to by p0. If p0 is null, a new random starting point is generated. If the
	// growth is successful, it returns a pointer to the subchain, or a NULL
	// pointer otherwise. It also sets the status variable to 1 for a fully
	// grown subchain, 0 otherwise.

	int		grown, loop;
	real		v[3];
	particleMD	p1, *pNew=0;

	// return if there is no monomer to add
	if (n==0) {
		*status = 1;
		return 0;
	}

	// add a monomer in the chain
	grown = 0;
	loop  = GROWLOOP_MAX;
	while (!grown && loop--) {

		// new monomer location
		if (p0) {
			RandomVector3D (v);
			p1.rx = p0->rx + v[x_]*sim->r0Fene/1.5;
			p1.ry = p0->ry + v[y_]*sim->r0Fene/1.5;
			// To allow 2d and 3d operation
			if (sim->box[z_] != 0.0 ){
				p1.rz = p0->rz + v[z_]*sim->r0Fene/1.5;
			}
			else {
				p1.rz = 0.0;
			}
			pNew = AtomInsert (sim, type, layout, &p1, CHECK, CHECK);
		}
		else {
			pNew = AtomInsert (sim, type, layout, 0, CHECK, CHECK);
		}

		// continue growing (recursively), and remove candidate if stunted growth
		if (pNew) {
			pNew->prev = p0;
			pNew->next = GrowLinearChain (sim, type, layout, n-1, pNew, &grown);
			if (!grown) AtomRemove (sim, type, PICK_POINTER, pNew);
		}
	}

	// update growth status
	*status = grown;

	// growth failure
	if (!grown) return NULL;

	// growth success
	return pNew;
}

//================================================================================
particleMD *GrowRodChain (simptr sim, int type, int layout, int n, particleMD *p0, int *status, int dir)
//================================================================================
{
	// RECURSIVE function to grow a n-monomer subchain from the particle pointed
	// to by p0. If p0 is null, a new random starting point is generated. If the
	// growth is successful, it returns a pointer to the subchain, or a NULL
	// pointer otherwise. It also sets the status variable to 1 for a fully
	// grown subchain, 0 otherwise.

	int		grown, loop;
	particleMD	p1, *pNew=0;

	// return if there is no monomer to add
	if (n==0) {
		*status = 1;
		return 0;
	}

	// add a monomer in the chain
	grown = 0;
	loop  = GROWLOOP_MAX;
	while (!grown && loop--) {

		// new monomer location
		if (p0) {
			p1.rx = p0->rx + (1-dir)*sim->r0Fene/1.5;
			p1.ry = p0->ry + (dir)*sim->r0Fene/1.5;
			p1.rz = p0->rz;
			pNew = AtomInsert (sim, type, layout, &p1, CHECK, CHECK);
		}
		else {
			pNew = AtomInsert (sim, type, layout, 0, CHECK, CHECK);
			// Force it to be at a give position rather than the random position it was inserted at
			pNew->rx = sim->box[x_]*0.5 - (1-dir)*n/2;
			pNew->ry = sim->box[y_]*0.5	- (dir)*n/2;
			pNew->rz = sim->box[z_]*0.5;
			pNew->wx = sim->box[x_]*0.5 - (1-dir)*n/2;
			pNew->wy = sim->box[y_]*0.5 - (dir)*n/2;
			pNew->wz = sim->box[z_]*0.5;
			pNew->x0 = sim->box[x_]*0.5 - (1-dir)*n/2;
			pNew->y0 = sim->box[y_]*0.5 - (dir)*n/2;
			pNew->z0 = sim->box[z_]*0.5;
		}

		// continue growing (recursively), and remove candidate if stunted growth
		if (pNew) {
			pNew->prev = p0;
			pNew->next = GrowRodChain (sim, type, layout, n-1, pNew, &grown, dir);
			if (!grown) AtomRemove (sim, type, PICK_POINTER, pNew);
		}
	}

	// update growth status
	*status = grown;

	// growth failure
	if (!grown) return NULL;

	// growth success
	return pNew;
}

//================================================================================
particleMD *GrowUChain (simptr sim, int type, int layout, int n, particleMD *p0, int *status)
//================================================================================
{
	// RECURSIVE function to grow a n-monomer subchain from the particle pointed
	// to by p0. If p0 is null, the starting point is calculated from the total
	// number of monomers. If the growth is successful, it returns a pointer to
	// the subchain, or a NULL pointer otherwise. It also sets the status
	// variable to 1 for a fully grown subchain, 0 otherwise.
	int		grown, loop;
	particleMD	p1, *pNew=0;

	// return if there is no monomer to add
	if (n==0) {
		*status = 1;
		return 0;
	}

	// add a monomer in the chain
	grown = 0;
	loop  = GROWLOOP_MAX;
	//even length monomer
	if ((sim->polyN[0] % 2) == 0) { // if even length
		while (!grown && loop--) {

			// new monomer location
			if (p0) {
				if (n > (sim->polyN[0]/2)) { // before turn
					p1.rx = p0->rx + sim->r0Fene/1.5;
					p1.ry = p0->ry;
					p1.rz = p0->rz;
				}
				else if (n == (sim->polyN[0]/2)) { // turn
					p1.rx = p0->rx;
					p1.ry = p0->ry - sim->r0Fene/1.5;
					p1.rz = p0->rz;
				} else { // after turn
					p1.rx = p0->rx - sim->r0Fene/1.5;
					p1.ry = p0->ry;
					p1.rz = p0->rz;
				}

				pNew = AtomInsert (sim, type, layout, &p1, CHECK, CHECK);
			}
			else {
				pNew = AtomInsert (sim, type, layout, 0, CHECK, CHECK);
				// Force it to be at a give position rather than the random position it was inserted at
				pNew->rx = sim->box[x_]*0.5 - n*0.25*sim->r0Fene/1.5; // n can be used as this is always the first iteration with highest n
				pNew->ry = sim->box[y_]*0.5 + sim->r0Fene/3.0;
				pNew->rz = sim->box[z_]*0.5;
				pNew->wx = sim->box[x_]*0.5 - n*0.25*sim->r0Fene/1.5;
				pNew->wy = sim->box[y_]*0.5 + sim->r0Fene/3.0;
				pNew->wz = sim->box[z_]*0.5;
				pNew->x0 = sim->box[x_]*0.5 - n*0.25*sim->r0Fene/1.5;
				pNew->y0 = sim->box[y_]*0.5 + sim->r0Fene/3.0;
				pNew->z0 = sim->box[z_]*0.5;
			}

			// continue growing (recursively), and remove candidate if stunted growth
			if (pNew) {
				pNew->prev = p0;
				pNew->next = GrowUChain (sim, type, layout, n-1, pNew, &grown);
				if (!grown) AtomRemove (sim, type, PICK_POINTER, pNew);
			}
		}
	}
	else { // if odd length
		while (!grown && loop--) {

			// new monomer location
			if (p0) {
				if (n > ((sim->polyN[0] + 1)/2)) { // before turn
					printf("before\n");
					p1.rx = p0->rx + sim->r0Fene/1.5;
					p1.ry = p0->ry;
					p1.rz = p0->rz;
				}
				else if (n == ((sim->polyN[0] + 1)/2)) { // turn1
					printf("turn1\n");
					p1.rx = p0->rx + sim->r0Fene*sqrt(3.0)/3.0;
					p1.ry = p0->ry - sim->r0Fene/3.0;
					p1.rz = p0->rz;
				}
				else if (n == (((sim->polyN[0] + 1)/2) - 1)) { // turn2
					printf("turn2\n");
					p1.rx = p0->rx - sim->r0Fene*sqrt(3.0)/3.0;
					p1.ry = p0->ry - sim->r0Fene/3.0;
					p1.rz = p0->rz;
				} else { // after turn
					printf("after\n");
					p1.rx = p0->rx - sim->r0Fene/1.5;
					p1.ry = p0->ry;
					p1.rz = p0->rz;
				}


				pNew = AtomInsert (sim, type, layout, &p1, CHECK, CHECK);
			}
			else {
				pNew = AtomInsert (sim, type, layout, 0, CHECK, CHECK);
				// Force it to be at a give position rather than the random position it was inserted at
				pNew->rx = sim->box[x_]*0.5 - (n+1)*0.25*sim->r0Fene/1.5; // n can be used as this is always the first iteration with highest n
				pNew->ry = sim->box[y_]*0.5 + sim->r0Fene/3.0;
				pNew->rz = sim->box[z_]*0.5;
				pNew->wx = sim->box[x_]*0.5 - (n+1)*0.25*sim->r0Fene/1.5;
				pNew->wy = sim->box[y_]*0.5 + sim->r0Fene/3.0;
				pNew->wz = sim->box[z_]*0.5;
				pNew->x0 = sim->box[x_]*0.5 - (n+1)*0.25*sim->r0Fene/1.5;
				pNew->y0 = sim->box[y_]*0.5 + sim->r0Fene/3.0;
				pNew->z0 = sim->box[z_]*0.5;
			}

			// continue growing (recursively), and remove candidate if stunted growth
			if (pNew) {
				pNew->prev = p0;
				pNew->next = GrowUChain (sim, type, layout, n-1, pNew, &grown);
				if (!grown) AtomRemove (sim, type, PICK_POINTER, pNew);
			}
		}
	}

	// update growth status
	*status = grown;

	// growth failure
	if (!grown) return NULL;

	// growth success
	return pNew;
}

//================================================================================
particleMD *AtomInsert (simptr sim, int type, int layout, particleMD *p, int chkolap, int chkptrs)
//================================================================================
{
	// Inserts an atom of the given type inside the atom list. This function
	// MUST be called before any other lists are created (e.g., neighbor lists,
	// anchor lists, charge lists etc.), because it only adjusts the main atom
	// list. Returns a pointer to the new atom upon success, and a NULL pointer
	// otherwise.

	int			picked, loop;
	int			i, n, n0, nType, nAtom;
	real		x=0, y=0, z=0;
	real		*box;
	particleMD 	*atom, *newAtom, *ptr;

	// local sim variables
	atom	= sim->atom.items;
	nAtom   = sim->atom.n;
	box	 	= sim->box;
	nType   = sim->nType[type];
	n0 		= sim->typeStart[type];

	// specific location (passed via non-null p)
	if (p) {
		x = p->rx;
		y = p->ry;
		z = p->rz;
		// enforce layout rule
		if (!LayoutRule(sim, layout, x, y, z)) return NULL;
		// reject overlap
		if (chkolap) if (AtomCheckOverlap(sim, type, x, y, z)) return NULL;
	}

	// random location (p is null)
	else {
		picked = 0;
		loop   = LOOP_MAX;
		while (!picked && loop--) {
			// pick a random location in box
			x = box[x_]*RandomReal();
			y = box[y_]*RandomReal();
			z = box[z_]*RandomReal();
			// enforce layout rule
			if (!LayoutRule(sim, layout, x, y, z)) continue;
			// reject overlap
			if (chkolap) if (AtomCheckOverlap(sim, type, x, y, z)) continue;
			picked=1;
		}
		if (!picked) error (ELOOPMAX);
	}

	// new atom
	newAtom = (particleMD *) mycalloc (1, sizeof(particleMD));
	newAtom->rx = newAtom->wx = newAtom->x0 = x;
	newAtom->ry = newAtom->wy = newAtom->y0 = y;
	newAtom->rz = newAtom->wz = newAtom->z0 = z;
	newAtom->type = type;
	newAtom->mass = sim->atomMass[type];
	newAtom->kspring = sim->atomSpring[type];
	newAtom->anchor = sim->atomAnchor[type];

	// calculate insertion index
	n = n0+nType;

	// insert at end of the type sublist
	if (n>0) {
		ptr = &(atom[n]);
		memmove ((ptr+1), ptr, (nAtom-n)*sizeof(particleMD));
		memmove (ptr, newAtom, sizeof(particleMD));
	}
	else {
		n = nAtom;
		sim->typeStart[type] = n;
		ptr = &(atom[n]);
		memmove (ptr, newAtom, sizeof(particleMD));
	}

	// adjust particle pointers
	if (chkptrs) OffsetParticlePointers (sim, n, 1);

	// adjust counts
	sim->nType[type]++;
	sim->atom.n++;
	nAtom++;

	// allocate additional memory for atoms if needed
	if (sim->atom.n >= sim->atom.max) GrowListAtom (&sim->atom);

	// adjust other types starting positions
	for (i=0; i<TYPE_COUNT; i++) {
		if (i==type) continue;
		if (sim->typeStart[i] >= n) sim->typeStart[i]++;
	}

	// free temporary space
	free (newAtom);

	return ptr;
}


//================================================================================
int AtomRemove (simptr sim, int type, int n, particleMD *p)
//================================================================================
{
	// Removes from the atom list the the atom pointed to by p (if p is not
	// null), or the n-th atom of a given type (n>0), a random atom of the given
	// type. This function MUST be called before any other lists are created
	// (e.g., neighbor lists, anchor lists, charge lists etc.), because it only
	// adjusts the main atom list. If p is not NULL, the particle pointed to by
	// p will be revomed. If the supplied index is negative, the function will
	// remove a random atom of the given type. Returns 1 if the atom was
	// indeed removed, 0 otherwise.

	int			picked, loop;
	int			i, n0, nType;
	particleMD 	*atom, *ptr;

	// get type from particle pointer
	if (p) type = p->type;

	// local sim variables
	atom	= sim->atom.items;
	nType	= sim->nType[type];
	n0	 	= sim->typeStart[type];

	// remove specific atom
	if (p) {
		// calculate deletion index and pointer
		n   = (int) (p - atom);
		ptr = p;
	}

	// random atom
	if (!p) {
		// return if there is no atom of this type
		if (nType==0) return 0;
		// specific index n
		if (n>=0) {
			// reject invalid index
			if (n >= nType) return 0;
			if (atom[n0+n].type != type) return 0;
		}
		// pick random index n
		else {
			picked = 0;
			loop   = LOOP_MAX;
			while (!picked && loop--) {
				// pick random atom index
				n = (int) (nType*RandomReal());
				// reject invalid index
				if (n == nType) continue;
				// we got a good one
				picked = 1;
			}
			if (!picked) error (ELOOPMAX);
		}
		// calculate deletion index and pointer
		n   = n0+n;
		ptr = &(atom[n]);
	}

	// adjust particle pointers (before deletion)
	OffsetParticlePointers (sim, n, -1);

	// remove atom
	memmove (ptr, ptr+1, (sim->atom.n-n-1)*sizeof(particleMD));

	// adjust counts
	sim->nType[type]--;
	if (sim->nType[type]==0) sim->typeStart[type] = 0;
	sim->atom.n--;

	// adjust other types starting positions
	for (i=0; i<TYPE_COUNT; i++) {
		if (i==type) continue;
		if (sim->typeStart[i] > n) sim->typeStart[i]--;
	}

	return 1;
}


//================================================================================
void OffsetParticlePointers (simptr sim, int n, int offset)
//================================================================================
{
	// This function updates the particle pointers when atoms are inserted or
	// removed from the atom list. All the pointers to locations after the
	// insertion or deletion point must be adjusted. Offset specifies offset in
	// number of particles; it is normally +1 for single atom insertion and -1
	// for single atom deletion.  Warning: pointers to a deleted atom are set to
	// NULL, i.e., bonds are lost.

	int		i, j, nAtom, nPolymer;
	particleMD	*atom, *p1, *p2;
	itemPoly	*polymer;

	// local sim variables
	atom	 = sim->atom.items;
	nAtom	 = sim->atom.n;
	polymer  = sim->polymer.items;
	nPolymer = sim->polymer.n;

	// return if there is no offset
	if (offset==0) return;

	// particle involved (insertion or deletion)
	p1 = atom+n;
	p2 = p1 + abs(offset) - 1;

	// adjust fene link pointers
	for (i=0; i<nAtom; i++) {
		if (atom[i].prev) OffsetPointer (&(atom[i].prev), p1, p2, offset);
		if (atom[i].next) OffsetPointer (&(atom[i].next), p1, p2, offset);
	}

	// adjust all polymer list pointers
	for (i=0; i<nPolymer; i++) {
		if (polymer[i].p1) OffsetPointer (&(polymer[i].p1), p1, p2, offset);
	}

	// adjust polymer list
	for (i=0; i<nPolymer; i++) {
	if (!polymer[i].p1) {
		for (j=i; j<nPolymer-1; j++) {
			polymer[j] = polymer[j+1];
		}
		nPolymer--;
		sim->polymer.n--;
	}
	}
}


//================================================================================
void OffsetPointer (particleMD **p, particleMD *p1, particleMD *p2, int offset)
//================================================================================
{
	// Offset the particle pointer pointed to by ptr, given the bounds p1 and p2
	// of particles involved, and the specified offset.

	// removing
	if (offset < 0) {
		if 		(*p>=p1 && *p<=p2) 	*p  = NULL;
		else if (*p > p2) 	   		*p += offset;
	}
	// inserting
	else {
		if (*p >= p1) *p += offset;
	}
}


/// Implements conditions to keep atoms or not when inserting atoms in the system.
/// The geometry of the system is thus really determined by this function. To add
/// new geometries, simply add their layout rules here.
///
/// @param		sim a pointer to a simulation structure
/// @param		layout the layout rule for the particle to insert
/// @param		x,y,z the point of insertion
/// @return		The function returns 1 if the atom is accepted for insertion,
///				and 0 otherwise.

//================================================================================
int LayoutRule (simptr sim, int layout, real x, real y, real z)
//================================================================================
{
	real r2=0, rMax2=0;

	// apply layout rule depending on geometry and type
	switch (sim->geometry) {

		case GEOM_BULK:
			return 1;
			break;

		case GEOM_PLATES:
			if ( y>sim->box[y_] || y<0) return 0;
			return 1;
			break;


		case GEOM_CYLINDER:
			// calculate distance from capillary axis (x)
			r2 =  (y-sim->caprIn-1.0)*(y-sim->caprIn-1.0) + (z-sim->caprIn-1.0)*(z-sim->caprIn-1.0);
			// rule for each layout type
			switch (layout) {
				case LAYOUT_WALL:
					// keep only if in cylindrical shell
					if (r2 > sim->caprIn2 && r2 < sim->caprOut2) return 1;
					else {
						if (r2 <= sim->caprIn2) sim->nPore+=1;
						return 0;
					}
					break;
				case LAYOUT_FLUID:
					// keep only if inside capillary
					rMax2  = sim->caprIn-0.5;
					rMax2 *= rMax2;
					if (r2 < rMax2) return 1;
					else return 0;
					break;
				//Tyler added. Copied from above
				case LAYOUT_RODX:
					// keep only if inside capillary
					rMax2  = sim->caprIn-0.5;
					rMax2 *= rMax2;
					if (r2 < rMax2) return 1;
					else return 0;
					break;
				//Zahra added. Copied from above
				case LAYOUT_RODY:
					// keep only if inside capillary
					rMax2  = sim->caprIn-0.5;
					rMax2 *= rMax2;
					if (r2 < rMax2) return 1;
					else return 0;
					break;
				// added by Karolina - identical to above
				case LAYOUT_U:
					// keep only if inside capillary
					rMax2  = sim->caprIn-0.5;
					rMax2 *= rMax2;
					if (r2 < rMax2) return 1;
					else return 0;
					break;
				case LAYOUT_CYLINDER:
					// keep only if inside capillary
					rMax2  = sim->caprIn;
					rMax2 *= rMax2;
					if (r2 < rMax2) return 1;
					else return 0;
					break;
				default:
				  printf("%d\n",layout);
					error (ELAYOUT);
			}
			break;

		case GEOM_CAPILLARY:
			// calculate distance from capillary axis (x)
			r2 =  y*y + z*z;
			// rule for each layout type
			switch (layout) {
				case LAYOUT_WALL:
					// keep only if in cylindrical shell
					if (r2 > sim->caprIn2 && r2 < sim->caprOut2) return 1;
					else {
						if (r2 <= sim->caprIn2) sim->nPore+=1;
						return 0;
					}
					break;
				case LAYOUT_FLUID:
					// keep only if inside capillary
					rMax2  = sim->caprIn-0.5;
					rMax2 *= rMax2;
					if (r2 < rMax2) return 1;
					else return 0;
					break;
				//Tyler added. Copied from above
				case LAYOUT_RODX:
					// keep only if inside capillary
					rMax2  = sim->caprIn-0.5;
					rMax2 *= rMax2;
					if (r2 < rMax2) return 1;
					else return 0;
					break;
				//Zahra added. Copied from above
				case LAYOUT_RODY:
					// keep only if inside capillary
					rMax2  = sim->caprIn-0.5;
					rMax2 *= rMax2;
					if (r2 < rMax2) return 1;
					else return 0;
					break;
				// added by Karolina - same as above
				case LAYOUT_U:
					// keep only if inside capillary
					rMax2  = sim->caprIn-0.5;
					rMax2 *= rMax2;
					if (r2 < rMax2) return 1;
					else return 0;
					break;
				default:
					error (ELAYOUT);
			}
			break;

		default:
			error (EGEOMETRY);
	}

	// don't keep by default
	return 0;
}


//================================================================================
int AtomCheckOverlap (simptr sim, int type, real x, real y, real z)
//================================================================================
{
	// Checks to see if an atom inserted at position x, y, z is within a
	// predefined distance of another atom in the system. Returns 1 if it is, 0
	// otherwise. This is a heavy function because all atom pairs are
	// considered.

	int			i;
	int			nAtom;
	real		dx, dy, dz, r2;
	real		olap1, olap2;
	particleMD 	*atom, *p;
	real		*box, *boxHalf;

	// local sim variables
	atom  = sim->atom.items;
	nAtom = sim->atom.n;
	olap1 = sim->overlapMin*sim->sigma_lj;
	olap2 = sim->overlapMinMonomer*sim->sigma_lj;
	box   = sim->box;
	boxHalf = sim->boxHalf;

	// square overlap minima
	olap1 *= olap1;
	olap2 *= olap2;


	// reject if there is overlap with other atoms
	for (i=0; i<nAtom; i++) {
		p = atom+i;

		// calculate position difference
		dx = p->rx - x;
		dy = p->ry - y;
		dz = p->rz - z;
		// Hard code the periodic BCs in
		if (sim->pbcond & PBC_COND_x) {
			if (dx >= boxHalf[x_]) dx -= box[x_];
			else if (dx < -boxHalf[x_]) dx += box[x_];
		}

		if (sim->pbcond & PBC_COND_y) {
			if (dy >= boxHalf[y_]) dy -= box[y_];
			else if	(dy < -boxHalf[y_]) dy += box[y_];
		}

		if (sim->pbcond & PBC_COND_z) {
			if (dz >= boxHalf[z_]) dz -= box[z_];
			else if	(dz < -boxHalf[z_]) dz += box[z_];
		}
//		ApplyPBC (sim, &dx, &dy, &dz);

		// olap1: overlap between all atoms
		r2 = dx*dx + dy*dy + dz*dz;
		if (r2 < olap1) {
			return 1;
		}

		// olap2: overlap for monomers (bigger to avoid tight knots)
		if (type == TYPE_MONOMER && p->type == TYPE_MONOMER) {
			if (r2 < olap2) {
				return 1;
			}
		}
		if (r2<0.2025) printf ("weird overlap, maybe extra parts floating? %E %d %d %E\n", r2, type, p->type, olap2);
	}
	return 0;
}
