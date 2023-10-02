///
///@file
///@brief 
///
/// Main loop that relies on subroutines to perform MPCD.

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* *********** SIMULATES AN MPCD GAS ********* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/*
By Tyler Shendruk's Research Group
	For the Active and Intelligent Matter Research Group
	At the University of Edinburgh

	Started late October 2008 at the University of Ottawa
	Continued at the University of Oxford
	Currently maintained at the University of Edinburgh
*/

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ LIBRARY FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
# include<stdio.h>
# include<math.h>
# include<stdlib.h>
# include<time.h>
# include<string.h>
#define _GNU_SOURCE // required for fenv
# include<fenv.h>
# include<signal.h>
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ DECLARE FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
# include "headers/definitions.h"
# include "headers/globals.h"
# include "headers/SRDclss.h"

# include "headers/mtools.h"
# include "headers/ctools.h"
# include "headers/rand.h"
# include "headers/read.h"
# include "headers/pout.h"
# include "headers/init.h"
# include "headers/mpc.h"
# include "headers/bc.h"
# include "headers/therm.h"
# include "headers/lc.h"
# include "headers/swimmers.h"

# include "../md/mdtypes.h"
# include "../md/mdsrd.h"
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* *************** MAIN THREAD ************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
///@brief
///
/// Call other functions to implement MPCD.
///
/// Declares all the program variable, reads the input file, initializes the system, perform a warm-up loop, then the main loop with checkpoints, before outputting and closing down.
///
/// @param argc Number of c argument variables.
/// @param argv[] Value of each c argument variables.
int main(int argc, char* argv[]) {
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	/* ************ PROGRAM VARIABLES *********** */
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	int i=0, j=0;					//Counting variables
	simptr simMD=NULL;				//The md simulation is a pointer
	cell ***CL;						//The array of all the cell lists
	particleMPC *SRDparticles;		//The array of particles
	spec *SPECIES;					//SPECIES is an array of the population of each species
	bc *WALL;						//The boundary conditions
	int starttime=0,runtime=0,warmtime=0;	//Temperal loop counter --- iterations NOT sim time
	//Input
	inputList inputVar;
	kinTheory theory;				//Theoretical values based on input
	//Timer variables
	time_t to,tf,lastCheckpoint;
	clock_t co,cf;
	//Simulation variables
	double KBTNOW=0.0;			//Current un-thermostated temperature
	double AVVEL=0.0;				//The average speed
	double AVS=0.0,S4=0.0,stdN=0.0;   //The average of the scalar order parameter, director and fourth moment
	double avDIR[_3D],AVV[_3D],AVNOW[_3D];  //The past and current average flow velocities
    zerovec_v(3, _3D, avDIR, AVV, AVNOW); // initialise to zero
	//Input/Output
	int CHCKPNTrcvr = 0;			//Flag for simulation from recovery of checkpoint
	char ip[500],op[500];			//Path to input and output
	int inMode = 0;					//Input mode: 0 - JSON, 1 - Legacy .inp
	outputFlagsList outFlags;		//Flags for what is outputted
	outputFilesList outFiles;		//List of output files
	specSwimmer specS;				//Swimmer's species
	swimmer *swimmers;				//Swimmers
	int WMD = 0;					// if MD works during warmup of mpcd

    /* ****************************************** */
    /* ****************************************** */
    /* ****************************************** */
    /* ************ Prepare Signals ************* */
    /* ****************************************** */
    /* ****************************************** */
    /* ****************************************** */
    #ifdef FPE
    #ifdef __linux__
    feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);
    #else
    printf("Floating point exception handling is only supported on Linux.\n");
    #endif
    #endif

	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	/* ************ BEGIN SIMULATION ************ */
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG > DBGRUN ) printf( "\nBegin SRD algorithm\n" );
	#endif
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	/* ************* READ ARGUMENTS ************* */
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG > DBGRUN ) printf("Read arguments and input files\n");
	#endif
	readarg( argc,argv,ip,op,&inMode );
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	/* ************ READ INPUT FILES ************ */
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	// read in depending on inMode
	if (inMode == 0){ // JSON input
		readJson( ip, &inputVar, &SPECIES, &SRDparticles, &CL, &MDmode, 
			&outFlags, &WALL, &specS, &swimmers);
	} else if (inMode == 1){ // Legacy .inp input
		readin( ip, &inputVar, &SPECIES, &SRDparticles, &CL, &MDmode );
		readbc( ip, &WALL );
		readpc( ip, &outFlags );
		readswimmers( ip, &specS, &swimmers );
	}	

	//Check if recovering checkpointed simulation
	if( inputVar.seed==-1 ) {
		CHCKPNTrcvr=1;
		#ifdef DBG
			if( DBUG >= DBGINIT ) printf( "Recovering checkpointed simulation\n" );
		#endif
		//Recovering checkpointed simulation
		readchckpnt( ip, &inputVar, &SPECIES, &SRDparticles, &CL, &MDmode, &WALL, &outFlags, &runtime, &warmtime, &theory, &AVVEL, &AVS, avDIR, &S4, &stdN, &KBTNOW, AVV, AVNOW, &specS, &swimmers );
	}
	#ifdef DBG
		if( DBUG > DBGRUN ) printf("Initialize simulation\n");
	#endif
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	/* ******** INITIALIZE OUTPUT FILES ********* */
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	initOutput( op, &outFlags, &outFiles, inputVar, SPECIES, WALL );
	if( outFlags.SYNOUT == OUT ) fprintf( outFiles.fsynopsis,"\nOutput files initialized.\n" );

	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	/* ************ INTIALIZE SYSTEM ************ */
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	// Check the input given by the user
	checkSim( outFiles.fsynopsis,outFlags.SYNOUT,inputVar,SPECIES,WALL,specS );

	if( CHCKPNTrcvr ) {
		//Simulation recovered from checkpoint
		#ifdef DBG
			if( DBUG >= DBGINIT ) printf("Simulation recovery initialization\n");
		#endif
		if( outFlags.SYNOUT == OUT ) fprintf( outFiles.fsynopsis,"Begin recovery initialization\n" );
		initializeRecovery( CL, SRDparticles, SPECIES, specS, inputVar.RTECH, inputVar.LC, MDmode, outFlags.SYNOUT, outFiles.fsynopsis );
	}
	else {
		#ifdef DBG
			if( DBUG >= DBGINIT ) printf("Initialize simulation variables\n");
		#endif
		if( outFlags.SYNOUT == OUT ) fprintf( outFiles.fsynopsis,"Initialize simulation variables\n" );
		// Initialize MD
		if( MDmode != noMD ) {
			#ifdef DBG
				if( DBUG >= DBGINIT ) printf( "\tLaunch MD simulation\n" );
			#endif
			simMD = launchMD(argc,argv);
			if(outFlags.SYNOUT == OUT) fprintf(outFiles.fsynopsis,"\nMD integrator launched.\n" );
		}
		//Normal initialization
		initializeSIM( CL, SRDparticles, SPECIES, WALL, simMD, &specS, swimmers, argc, argv, &inputVar, &to, &co, &runtime, &warmtime, &AVVEL, &theory, &KBTNOW, &AVS, &S4, &stdN, AVNOW, AVV, avDIR, outFlags, MDmode, outFiles.fsynopsis, ip );
	}
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	/* ***************** INPUT ****************** */
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	stateinput( inputVar, SPECIES, WALL, specS, outFlags, theory, outFiles.fsynopsis );
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	/* *************** STATE INPUT ************** */
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	listinput( inputVar, AVVEL, SPECIES, theory );
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	/* **************** MAIN BODY *************** */
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
    lastCheckpoint = time(NULL); // set NOW to last checkpoint time

	/* ****************************************** */
	/* *************** WARMUP LOOP ************** */
	/* ****************************************** */
	if( inputVar.warmupSteps>0 && !CHCKPNTrcvr ) {
		#ifdef DBG
			if( DBUG > DBGRUN ) printf( "Begin warmup loop\n" );
		#endif
		if(outFlags.SYNOUT == OUT) fprintf( outFiles.fsynopsis,"\nBegin warmup loop.\n" );
		// This is the main loop of the SRD program. The temporal loop. 
		starttime=warmtime;
		if(simMD->warmupMD){
			WMD = 1;
		}
		for( warmtime=starttime; warmtime<=inputVar.warmupSteps; warmtime++ ) {
			/* ****************************************** */
			/* ***************** UPDATE ***************** */
			/* ****************************************** */
			timestep( CL, SRDparticles, SPECIES, WALL, simMD, &specS, swimmers, AVNOW, AVV, avDIR, inputVar, &KBTNOW, &AVS, warmtime, MDmode, outFlags, outFiles ,WMD);
			/* ****************************************** */
			/* *************** CHECKPOINT *************** */
			/* ****************************************** */
			if( outFlags.CHCKPNT>=OUT && warmtime%outFlags.CHCKPNT==0 ) {
                runCheckpoint( op, &lastCheckpoint, outFiles.fchckpnt, inputVar, SPECIES, SRDparticles, MDmode, WALL, outFlags, runtime, warmtime, AVVEL, AVS, avDIR, S4, stdN, KBTNOW, AVV, AVNOW, theory,specS, swimmers );
			}
		}
		inputVar.warmupSteps=0;
		#ifdef DBG
			if( DBUG > DBGRUN ) printf( "Warmup loop complete\n" );
		#endif
		if(outFlags.SYNOUT == OUT) fprintf( outFiles.fsynopsis,"\nWarmup loop complete\n" );
	}
	/* ****************************************** */
	/* ************** TEMPORAL LOOP ************* */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG > DBGRUN ) printf( "Begin temporal loop\n" );
	#endif
	if(outFlags.SYNOUT == OUT) fprintf( outFiles.fsynopsis,"\nBegin temporal loop.\n" );
	// This is the main loop of the SRD program. The temporal loop.
	starttime=runtime;
	WMD = 1;
	for( runtime=starttime; runtime<=inputVar.simSteps; runtime++ ) {
		/* ****************************************** */
		/* ***************** UPDATE ***************** */
		/* ****************************************** */
		timestep( CL, SRDparticles, SPECIES, WALL, simMD, &specS, swimmers, AVNOW, AVV, avDIR, inputVar, &KBTNOW, &AVS, runtime, MDmode, outFlags, outFiles,WMD );
		/* ****************************************** */
		/* ***************** OUTPUT ***************** */
		/* ****************************************** */
		if( writeOutput(runtime,outFlags,inputVar.RFRAME,inputVar.zeroNetMom) ) {
			outputResults( CL, SRDparticles, SPECIES, WALL, simMD, specS, swimmers, AVNOW, AVV, avDIR, runtime, inputVar, AVVEL, KBTNOW, &AVS, &S4, &stdN, MDmode, outFlags, outFiles );
		}
		//The histogram stuff takes a lot of memory so only make it IF necessary
		if( writeHistograms(runtime,outFlags) ) outputHist( CL, runtime, inputVar, outFlags, outFiles );
		/* ****************************************** */
		/* *************** CHECKPOINT *************** */
		/* ****************************************** */
		if( outFlags.CHCKPNT>=OUT && runtime%outFlags.CHCKPNT==0 ) {
            runCheckpoint( op, &lastCheckpoint, outFiles.fchckpnt, inputVar, SPECIES, SRDparticles, MDmode, WALL, outFlags, runtime, warmtime, AVVEL, AVS, avDIR, S4, stdN, KBTNOW, AVV, AVNOW, theory,specS, swimmers );
        }
	}
	#ifdef DBG
		if( DBUG > DBGRUN ) printf( "Temporal loop complete\n" );
	#endif
	if( outFlags.SYNOUT == OUT ) fprintf( outFiles.fsynopsis,"\nTemporal loop complete\n" );
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	/* ******************* END ****************** */
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	tf = time( NULL );
	cf = clock( );
	#ifdef DBG
		if( DBUG >= DBGINIT ) {
            float cpuTime = (float) (cf-co)/CLOCKS_PER_SEC;
            printf( "Wall compuation time: %e sec\nCPU compuation time:  %e CPUsec\n",(float)(tf-to), cpuTime );

            // compute particle updates per second (PUPS)
            int mpcUpdates = (inputVar.simSteps+inputVar.warmupSteps) * GPOP; // total # of mpc particle updates
            int mpcPUPS = mpcUpdates / cpuTime; // mpc particle updates per second

            // convert PUPS to KPUPS or MPUPs if necessary
            if (mpcPUPS > 1000000) {
                printf("\tMPCD Particle Updates Per Second: %.1f MPUPS\n\n", (float) mpcPUPS/1000000);
            } else if (mpcPUPS > 1000) {
                printf("\tMPCD Particle Updates Per Second: %.1f KPUPS\n\n", (float) mpcPUPS/1000);
            } else {
                printf("\tMPCD Particle Updates Per Second: %d PUPS\n\n", mpcPUPS);
            }
            //NOTE: This doesn't take MD into account
        }
	#endif
	if( outFlags.SYNOUT ) {
		fprintf( outFiles.fsynopsis, "\nWall compuation time: %e sec\nCPU compuation time:  %e CPUsec\n\n",(float)(tf-to),(float) (cf-co)/CLOCKS_PER_SEC );
		fclose( outFiles.fsynopsis );
	}
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "Closing Files.\n" );
	#endif
	closeOutputFiles( SPECIES, WALL, outFlags, outFiles );

	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "Freeing memory.\n" );
	#endif
	free( SRDparticles );
	free( SPECIES );
	free( WALL );
	free( swimmers );
	//To free the D3D array must free z-values, then free the y-pointers to those arrays then free the x-pointers to the array of y-pointers
	for( i=0; i<XYZ_P1[0]; i++ ) {
		for( j=0; j<XYZ_P1[1]; j++ ) free( CL[i][j] );
		free( CL[i] );
	}
	free( CL );
	#ifdef DBG
		if( DBUG > DBGRUN ) printf( "Simulation successful\n" );
		if( DBUG >= DBGWARN ) printf("\a");
	#endif
	return EXIT_SUCCESS;
}
