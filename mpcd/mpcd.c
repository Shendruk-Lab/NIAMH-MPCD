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
#if defined(_OPENMP)
#include <omp.h>
#endif
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
int main(int argc, char* argv[]) {
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	/* ************ PROGRAM VARIABLES *********** */
	/* ****************************************** */
	/* ****************************************** */
	/* ****************************************** */
	int i,j;						//Counting variables
	simptr simMD=NULL;				//The md simulation is a pointer
	cell ***CL;						//The array of all the cell lists
	particleMPC *SRDparticles;		//The array of particles
	spec *SPECIES;					//SPECIES is an array of the population of each species
	bc *WALL;						//The boundary conditions
	int starttime,runtime,warmtime;	//Temperal loop counter --- iterations NOT sim time
	//Input
	inputList inputVar;
	kinTheory theory;				//Theoretical values based on input
	//Timer variables
	time_t to,tf,lastCheckpoint;
	clock_t co,cf;
	//Simulation variables
	double KBTNOW;					//Current un-thermostated temperature
	double AVVEL;					//The average speed
	double AVS,avDIR[_3D],S4,stdN;	//The average of the scalar order parameter, director and fourth moment
	double AVV[_3D],AVNOW[_3D];		//The past and current average flow velocities
	//Input/Output
	int CHCKPNTrcvr = 0;			//Flag for simulation from recovery of checkpoint
	char ip[500],op[500];			//Path to input and output
	int inMode = 0;					//Input mode: 0 - JSON, 1 - Legacy .inp
	outputFlagsList outFlags;		//Flags for what is outputted
	outputFilesList outFiles;		//List of output files
	specSwimmer specS;				//Swimmer's species
	swimmer *swimmers;				//Swimmers

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
#ifdef _OPENMP
#ifdef DBG
    ///TODO: This should have a "debug level" check but I can't get it to cooperate
    printf("OpenMP is enabled with %d max threads\n", omp_get_max_threads());
#endif
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
		for( warmtime=starttime; warmtime<=inputVar.warmupSteps; warmtime++ ) {
			/* ****************************************** */
			/* ***************** UPDATE ***************** */
			/* ****************************************** */
			timestep( CL, SRDparticles, SPECIES, WALL, simMD, &specS, swimmers, AVNOW, AVV, avDIR, inputVar, &KBTNOW, &AVS, warmtime, MDmode, outFlags, outFiles );
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
	/* ************** TEMPERAL LOOP ************* */
	/* ****************************************** */
	#ifdef DBG
		if( DBUG > DBGRUN ) printf( "\nBegin temperal loop\n" );
	#endif
	if(outFlags.SYNOUT == OUT) fprintf( outFiles.fsynopsis,"\nBegin temperal loop.\n" );
	// This is the main loop of the SRD program. The temporal loop.
	starttime=runtime;
	for( runtime=starttime; runtime<=inputVar.simSteps; runtime++ ) {
		/* ****************************************** */
		/* ***************** UPDATE ***************** */
		/* ****************************************** */
		timestep( CL, SRDparticles, SPECIES, WALL, simMD, &specS, swimmers, AVNOW, AVV, avDIR, inputVar, &KBTNOW, &AVS, runtime, MDmode, outFlags, outFiles );
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
		if( DBUG > DBGRUN ) printf( "Temperal loop complete\n" );
	#endif
	if( outFlags.SYNOUT == OUT ) fprintf( outFiles.fsynopsis,"\nTemperal loop complete\n" );
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
		if( DBUG >= DBGINIT ) printf( "Wall compuation time: %e sec\nCPU compuation time:  %e CPUsec\n\n",(float)(tf
			-to),(float) (cf-co)/CLOCKS_PER_SEC );
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
