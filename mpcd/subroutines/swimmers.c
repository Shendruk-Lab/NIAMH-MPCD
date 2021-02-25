# include<math.h>
# include<stdio.h>
# include<stdlib.h>
# include<string.h>

# include "../headers/SRDclss.h"
# include "../headers/definitions.h"
# include "../headers/globals.h"
# include "../headers/swimmers.h"
# include "../headers/rand.h"
# include "../headers/mtools.h"
# include "../headers/init.h"
# include "../headers/pout.h"
# include "../headers/mdbc.h"
# include "../headers/mpc.h"
# include "../headers/bc.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ******* ROUTINES FOR READING FILES ******* */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
void readswimmers( char fpath[],specSwimmer *specS,swimmer **sw ) {
/*
   By reading in the addresses of the variables as
   pointers this function sets the values to what
   is in the input file swimmer.inp
*/
	FILE *finput;
	char STR[100];
	int MS;
	double MF;

	strcpy( STR,fpath );
	strcat( STR,"swimmer.inp" );
	finput = fopen( STR, "r" );
	if( !finput ) {					// file couldn't be opened
		printf( "Error:\tFile '%s' could not be opened.\n",STR );
		exit( 1 );
	}

	if(fscanf( finput,"%d %s",&MS,STR ));		//Read swimmer type
	else printf("Warning: Failed to read swimmer typer.\n");
	specS->TYPE = MS;
	if(fscanf( finput,"%d %s",&MS,STR ));			//Read swimmer population
	else printf("Warning: Failed to read swimmer populations.\n");
	NS = MS;
	if(fscanf( finput,"%d %s",&MS,STR ));			//Read position initialization
	else printf("Warning: Failed to read swimmer position initialization.\n");
	specS->QDIST = MS;
	if(fscanf( finput,"%d %s",&MS,STR ));			//Read orientation initialization
	else printf("Warning: Failed to read swimmer orientation initialization.\n");
	specS->ODIST = MS;
	if(fscanf( finput,"%d %s",&MS,STR ));			//Read monomer head mass
	else printf("Warning: Failed to read swimmer-monomer mass.\n");
	specS->headM = MS;
	if(fscanf( finput,"%d %s",&MS,STR ));			//Read monomer tail mass
	else printf("Warning: Failed to read swimmer-monomer mass.\n");
	specS->middM = MS;
	if(fscanf( finput,"%d %s",&MS,STR ));			//Read head monomer type
	else printf("Warning: Failed to read swimmer-head monomer type.\n");
	specS->HSPid = MS;
	if(fscanf( finput,"%d %s",&MS,STR ));			//Read tail monomer type
	else printf("Warning: Failed to read swimmer-tail monomer type.\n");
	specS->MSPid = MS;
	if(fscanf( finput,"%lf %s",&MF,STR ));			//Read swimming propulsion force
	else printf("Warning: Failed to read swimming propulsion force.\n");
	specS->FS = MF;
  if(fscanf( finput,"%lf %s",&MF,STR ));			//Read swimming dipole strength (i.e. length of dipole in multiples of head/middle separation)
	else printf("Warning: Failed to read swimming dipole strength.\n");
	specS->DS = MF;
	if(fscanf( finput,"%lf %s",&MF,STR ));			//Read swimming rotlet dipole torque
	else printf("Warning: Failed to read swimming dipole strength.\n");
	specS->TS = MF;

	if(fscanf( finput,"%lf %s",&MF,STR ));			//Read tumbling shrink size (i.e. what ro and sigma are scaled by for the tumbling)
  else printf("Warning: Failed to read tumbling shrink size coefficient.\n");
  specS->sizeShrink = MF;
	if(fscanf( finput,"%lf %s",&MF,STR ));			//Read tumbling shrink spring (i.e. what k is scaled by for tumbling)
	else printf("Warning: Failed to read tumbling shrink spring coefficient.\n");
	specS->springShrink = MF;

	if(fscanf( finput,"%lf %s",&MF,STR ));			//Read swimmers' spring constant
	else printf("Warning: Failed to read swimmer-FENE spring constant.\n");
	specS->k = MF;
	if(fscanf( finput,"%lf %s",&MF,STR ));			//Read swimmers' spring constant
	else printf("Warning: Failed to read swimmer-FENE ro.\n");
	specS->ro = MF;
	if(fscanf( finput,"%lf %s",&MF,STR ));		//Read swimmers' monomer size
	else printf("Warning: Failed to read swimmer-LJ size.\n");
	specS->sig = MF;
	if(fscanf( finput,"%lf %s",&MF,STR ));		//Read swimmers' energy
	else printf("Warning: Failed to read swimmer-LJ energy.\n");
	specS->eps = MF;

  if(fscanf( finput,"%lf %s",&MF,STR ));		//Read swimmers' average run time
	else printf("Warning: Failed to read swimmer average run time.\n");
	specS->runTime = MF;
	if(fscanf( finput,"%lf %s",&MF,STR ));		//Read swimmers' average tumble time
	else printf("Warning: Failed to read swimmer average tumble time.\n");
	specS->tumbleTime = MF;
	if(fscanf( finput,"%d %s",&MS,STR ));		//Read swimmers' shrink/extension time
	else printf("Warning: Failed to read swimmer shrink time.\n");
	specS->shrinkTime = MS;

	if(fscanf( finput,"%lf %s",&MF,STR ));		//Read swimmers' magnetic susceptibility
	else printf("Warning: Failed to read swimmer magnetic susceptibility.\n");
	specS->MAGMOM = MS;

	if(fscanf( finput,"%lf %s",&MF,STR ));			//Read swimming rotlet dipole torque
	else printf("Warning: Failed to read swimming dipole strength.\n");
	specS->fixDist = MF;

	//Allocate the memory for the swimmers
	(*sw) = (swimmer*) malloc( NS * sizeof( swimmer ) );

	fclose( finput );
}
void setswimmers( specSwimmer *SS,swimmer *swimmers,bc WALL[],int stepsMD,double dt ) {
/*
    This subroutine does the initializes the swimmers' coordinates.
*/
	int i,j,d,flag;
	double U[DIM],W;
  double shift[DIM];

	SS->iheadM=1.0/(double)SS->headM;
	SS->imiddM=1.0/(double)SS->middM;
	SS->iro=1.0/SS->ro;
	SS->isig=1.0/SS->sig;
	//Zero coordinates
	for( i=0; i<NS; i++ ) {
		for( d=0; d<_3D; d++ ) {
			(swimmers+i)->H.Q[d]=0.0;
			(swimmers+i)->H.V[d]=0.0;
			(swimmers+i)->H.A[d]=0.0;
			(swimmers+i)->M.Q[d]=0.0;
			(swimmers+i)->M.V[d]=0.0;
			(swimmers+i)->M.A[d]=0.0;
		}
		//Set head and middle
		(swimmers+i)->H.HorM=0;
		(swimmers+i)->M.HorM=1;
	}
	//Set coordinates
	if( SS->TYPE==DUMBBELL_FIXED ) {
		if(NS>1) {
			printf( "Error: Swimmer type DUMBBELL_FIXED is reserved for a single fixed dumbell. Please set swimmer number to 1.\n" );
			exit( 1 );
		}
		for( d=0; d<DIM; d++ ) swimmers->H.Q[d] = 0.5 * (float)XYZ[d];
		swimmers->H.Q[0] += 0.5*SS->ro;
		for( d=0; d<DIM; d++ ) swimmers->M.Q[d] = 0.5 * (float)XYZ[d];
		swimmers->M.Q[0] -= 0.5*SS->ro;
	}
	else {
	  for( i=0; i<NS; i++ ) {
	    flag=1;
	    while( flag ) {
	  		//Set coordinates of the swimmer's head
	  		if( SS->QDIST == PRF ) for( d=0; d<DIM; d++ ) (swimmers+i)->H.Q[d] = genrand_real( ) * XYZ[d];
	  		else if( SS->QDIST == CENTRED ) for( d=0; d<DIM; d++ ) (swimmers+i)->H.Q[d] = 0.5 * XYZ[d];
	  		else{
	  			printf( "Error: Swimmer placement distribution unacceptable.\n" );
	  			exit( 1 );
	  		}
				// orient() should produce the correct direction to do all ODIST options
	  		if( SS->ODIST == RANDORIENT ) {
	  			orient( U,SS->ODIST );
	  			for( d=0; d<DIM; d++ ) {
	  				(swimmers+i)->M.Q[d] = (swimmers+i)->H.Q[d] - 0.49*SS->ro*U[d];
	  				(swimmers+i)->H.Q[d] = (swimmers+i)->H.Q[d] + 0.49*SS->ro*U[d];
	  			}
	  		}
	      // Check if the placement is acceptable
				flag=0;
	      for( j=0; j<NBC; j++ ) {
	        // Head
	        shiftBC_swimmer( shift,&WALL[j],&(swimmers+i)->H );
	        W = calcW_swimmer( WALL[j],&(swimmers+i)->H );
	        if( W<=0 ) flag=1;
	        shiftbackBC( shift,&WALL[j] );
	        // Middle
	        shiftBC_swimmer( shift,&WALL[j],&(swimmers+i)->M );
	        W = calcW_swimmer( WALL[j],&(swimmers+i)->M );
	        if( W<=0 ) flag=1;
	        shiftbackBC( shift,&WALL[j] );
					// Don't allow PBCs
					for( d=0; d<DIM; d++ ) {
						if( (swimmers+i)->H.Q[d]>XYZ[d] ) flag=1;
						else if ( (swimmers+i)->H.Q[d]<0 ) flag=1;
						if( (swimmers+i)->M.Q[d]>XYZ[d] ) flag=1;
						else if ( (swimmers+i)->M.Q[d]<0 ) flag=1;
					}
	      }
				// If excluded volume then check other swimmers too
				if(SS->TYPE==DUMBBELL_EXVOL) for( j=0; j<i; j++ ) {
					for( d=0; d<_3D; d++ ) shift[d]=0.0;
					//Head-head
					for( d=0; d<DIM; d++ ) shift[d]=(swimmers+i)->H.Q[d]-(swimmers+j)->H.Q[d];
					if(length(shift,DIM)<SS->sig) flag=1;
					//Head-middle
					for( d=0; d<DIM; d++ ) shift[d]=(swimmers+i)->H.Q[d]-(swimmers+j)->M.Q[d];
					if(length(shift,DIM)<SS->sig) flag=1;
					//Middle-head
					for( d=0; d<DIM; d++ ) shift[d]=(swimmers+i)->M.Q[d]-(swimmers+j)->H.Q[d];
					if(length(shift,DIM)<SS->sig) flag=1;
					//Middle-middle
					for( d=0; d<DIM; d++ ) shift[d]=(swimmers+i)->M.Q[d]-(swimmers+j)->M.Q[d];
					if(length(shift,DIM)<SS->sig) flag=1;
				}
	    }
	  }
	}
  //Give each swimmer it's size and chain
  for( i=0; i<NS; i++ ) {
    (swimmers+i)->ro = SS->ro;
    (swimmers+i)->sig = SS->sig;
    (swimmers+i)->iro = 1.0/(swimmers+i)->ro;
    (swimmers+i)->isig = 1.0/(swimmers+i)->sig;
		(swimmers+i)->k = SS->k;
  }
  //Set run/tumble
  if( feq(SS->runTime,0.0) || feq(SS->tumbleTime,0.0) || SS->shrinkTime==0 ) {
    SS->runTime=0.0;
    SS->tumbleTime=0.0;
		SS->shrinkTime=0;
  }
  else for( i=0; i<NS; i++ ) {
    W=genrand_real();
    //Set the fraction running or tumbling according to their relative times
		//Running phase
    if( W<=SS->runTime/(SS->runTime+SS->tumbleTime) ) {
      (swimmers+i)->RT = RUNNING;   //Set swimmer in run phase
      (swimmers+i)->timeRND = (int) genrand_exp( SS->runTime );    //Generate random integer to run for
      (swimmers+i)->timeCNT = (int) (genrand_real()*((swimmers+i)->timeRND));    //Swimmer timer/counter set to zero
    }
		//Tumbling phase
    else{
      (swimmers+i)->RT = TUMBLING;   //Set swimmer in tumbling phase initially
      (swimmers+i)->timeRND = (int) genrand_exp( SS->tumbleTime );    //Generate random integer to tumble for
      (swimmers+i)->timeCNT = (int) (genrand_real()*((swimmers+i)->timeRND));    //Swimmer timer/counter set to zero
      //Shrink the swimmer
      if( fneq(SS->sizeShrink,1.0) || fneq(SS->springShrink,1.0) ) {
				#ifdef DBG
					if( DBUG == DBGRUNTUMBLE ) {
						printf( "\tInitial ro=%lf, sig=%lf, k=%lf\n",(swimmers+i)->ro,(swimmers+i)->sig,(swimmers+i)->k );
						printf( "\tShrink size coef=%lf, spring=%lf\n",SS->sizeShrink,SS->springShrink );
					}
				#endif
				(swimmers+i)->sig*=(SS->sizeShrink);
				(swimmers+i)->isig=1.0/((swimmers+i)->sig);
				(swimmers+i)->ro*=(SS->sizeShrink);
				(swimmers+i)->iro=1.0/((swimmers+i)->ro);
				(swimmers+i)->k*=(SS->springShrink);
				for( d=0; d<DIM; d++ ) shift[d] = (swimmers+i)->M.Q[d] - (swimmers+i)->H.Q[d];
				#ifdef DBG
					if( DBUG == DBGRUNTUMBLE ) {
						printf( "\tseparation: %lf, ",length(shift,DIM) );
						pvec(shift,DIM);
					}
				#endif
				for( d=0; d<DIM; d++ ) (swimmers+i)->M.Q[d] = (swimmers+i)->H.Q[d] + shift[d]*(SS->sizeShrink);
				#ifdef DBG
					if( DBUG == DBGRUNTUMBLE ) {
						printf( "\tFinal ro=%lf, sig=%lf, k=%lf\n",(swimmers+i)->ro,(swimmers+i)->sig,(swimmers+i)->k );
					}
				#endif
				#ifdef DBG
					if( DBUG == DBGRUNTUMBLE ) {
						for( d=0; d<DIM; d++ ) shift[d] = (swimmers+i)->M.Q[d] - (swimmers+i)->H.Q[d];
						printf( "\tseparation: %lf, ",length(shift,DIM) );
						pvec(shift,DIM);
					}
				#endif
      }
    }
  }
	// Calculate the swimmers' orientations
	for( i=0; i<NS; i++ ) {
		for( d=0; d<_3D; d++ ) (swimmers+i)->n0[d] = 0.0;
		swimmerOri( (swimmers+i)->n0,(swimmers+i) );
	}
}
void swimout( FILE *fout,swimmer swimmers[],double t ) {
/*
    Print swimmer positions data to file
*/
	int i;

	for( i=0; i<NS; i++ ) fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%d\n",t,swimmers[i].H.Q[0],swimmers[i].H.Q[1],swimmers[i].H.Q[2],swimmers[i].H.V[0],swimmers[i].H.V[1],swimmers[i].H.V[2],swimmers[i].M.Q[0],swimmers[i].M.Q[1],swimmers[i].M.Q[2],swimmers[i].M.V[0],swimmers[i].M.V[1],swimmers[i].M.V[2],swimmers[i].RT );
	fprintf( fout, "\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void swimoriout( FILE *fout,swimmer swimmers[],double t ) {
/*
    Print swimmer orientation data to file
*/
	int i,d;
	double n[_3D];

	for( d=0; d<_3D; d++ ) n[d]=0.0;
	for( i=0; i<NS; i++ ) {
		swimmerOri( n,(swimmers+i) );
		fprintf( fout, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%d\n",t,n[0],n[1],n[2],swimmers[i].RT );
	}
	fprintf( fout, "\n" );
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void runtumbleout( FILE *fout,swimmer *sw ) {
/*
    Print swimmer positions data to file
*/
	int i;
	double n[_3D];				//Swimmer orientation
	double dAng;

	for( i=0; i<_3D; i++ ) n[i]=0.0;
	for( i=0; i<NS; i++ ) {
		swimmerOri( n,sw );
		dAng=signedAngle( sw->n0,n,_3D );
		fprintf( fout,"%d\t%d\t%12.5e\n",sw->RT,sw->timeCNT,dAng );
	}
	#ifdef FFLSH
		fflush(fout);
	#endif
}
void swimmerPBC_dr(double *dr ) {
/*
    Apply period boundary conditions to swimmer monomer interactions
*/
  int d;
	double boxHalf[DIM];

	for( d=0; d<DIM; d++ ) boxHalf[d] = 0.5*XYZ[d];
	for( d=0; d<DIM; d++ ) {
		if( dr[d] >= boxHalf[d] ) dr[d] -= XYZ[d];
		else if( dr[d] < -boxHalf[d]) dr[d] += XYZ[d];
	}
}
void swimmerOri( double n[],swimmer *sw ) {
/*
   Orientation vector of a swimmer
*/
	int d;
	//Vector from middle to head
	for( d=0; d<DIM; d++ ) n[d] = sw->H.Q[d] - sw->M.Q[d];
	//Apply potential PBCs
	swimmerPBC_dr( n );
	//Normalize
	norm( n,DIM );
}
double swimmerWCA( double r,double eps ) {
/*
	Calculate the WCA force from the separation r SCALED by sigma=size
	To get the force vector this must be multiplied by vec(r)
*/
	double r2,r6,f;
	double rcut=1.122462048309373, fcap=1E3;

	r2 = 1./(r*r);
	r6 = r2*r2*r2;
	if( r<rcut ) {
		f = 48.0 * eps * r6 * r2 * (r6 - 0.5);
		if( f>fcap ) f=fcap;
		return f;
	}
	else return 0.0;
}
double swimmerFENE( double r,double k ) {
/*
	Calculate the FENE force from the separation r SCALED by ro=size
	To get the force vector this must be multiplied by vec(r)
	If the FENE chain is passed then there is a large Hookean "backup" force to pull them together
*/
	#ifdef DBG
			if( DBUG == DBGSWIMMER || DBUG == DBGSWIMMERDEETS ) if( r>=1.0 ) printf("Warning: dr=%e>1: FENE spring overstretched. Will use strong harmonic spring with 50k\n",r);
	#endif
	if( r<1.0 ) return k/(r*r-1.0);								//FENE
	// else return swimmerHookean( r,50.0*k );		//Linear hookean force
	else return swimmerSpring6( r,k );				//Non-linear spring force
}
double swimmerHookean( double r,double k ) {
/*
	Calculate the Hookean spring force from the separation r SCALED by ro=size
	To get the force vector this must be multiplied by vec(r)
*/
	return -k*r;
}
double swimmerSpring6( double r,double k ) {
/*
	Calculate the non-linear spring force from the separation r SCALED by ro=size
	To get the force vector this must be multiplied by vec(r)
*/
	double r2=r*r;
	return -k*r*r2*r2;
}
void swimmerVerlet_nonInteracting( specSwimmer SS,swimmer *s,double dt,int springType, cell ***CL) {
/*
	The leapfrog/Verlet-Velocity algorithm
*/
	int d;
	double a[DIM];

	//Verlet step 1
	//Update velocity using last time step's acceleration
	for( d=0; d<DIM; d++ ) s->H.V[d] += 0.5*dt*s->H.A[d];
	for( d=0; d<DIM; d++ ) s->M.V[d] += 0.5*dt*s->M.A[d];

	//Verlet step 2
	//Update position using the half-step velocity
	for( d=0; d<DIM; d++ ) s->H.Q[d] += dt*s->H.V[d];
	for( d=0; d<DIM; d++ ) s->M.Q[d] += dt*s->M.V[d];

	//Verlet step 3
	//Update the forces using the new positions
	//Calculate force on the head (and know force on middle is equal but opposite)
	smonoForce_sameSwimmer( a,SS,s,springType );
	//Save the acceleration
	for( d=0; d<DIM; d++ ) s->H.A[d] = a[d]*SS.iheadM;
	for( d=0; d<DIM; d++ ) s->M.A[d] = -1.0*a[d]*SS.imiddM;

	// //Same for the alingment potential 
	// smonoForce_AlignmentPot(a, CL,SS,(swimmers + i));
	// for( d=0; d<DIM; d++ ) (swimmers+i)->H.A[d] += a[d]*SS.iheadM;
	// for( d=0; d<DIM; d++ ) (swimmers+i)->M.A[d] -= a[d]*SS.imiddM;

	//Verlet step 4
	//Update velocity
	for( d=0; d<DIM; d++ ) s->H.V[d] += 0.5*dt*s->H.A[d];
	for( d=0; d<DIM; d++ ) s->M.V[d] += 0.5*dt*s->M.A[d];
}
void swimmerVerlet_all( specSwimmer SS,swimmer swimmers[],double dt,int springType, cell ***CL) {
/*
	The leapfrog/Verlet-Velocity algorithm
*/
	int i,j,d;
	double a[DIM]; //force

	//Verlet step 1
	//Update velocity using last time step's acceleration
	for( i=0; i<NS; i++ ) for( d=0; d<DIM; d++ ) {
		(swimmers+i)->H.V[d] += 0.5*dt*(swimmers+i)->H.A[d];
		(swimmers+i)->M.V[d] += 0.5*dt*(swimmers+i)->M.A[d];
	}

	//Verlet step 2
	//Update position using the half-step velocity (zero the acceleration so can update it in step 3)
	for( i=0; i<NS; i++ ) for( d=0; d<DIM; d++ ) {
		(swimmers+i)->H.Q[d] += dt*(swimmers+i)->H.V[d];
		(swimmers+i)->M.Q[d] += dt*(swimmers+i)->M.V[d];
		(swimmers+i)->H.A[d] = 0.0;
		(swimmers+i)->M.A[d] = 0.0;
	}

	//Verlet step 3
	//Update the forces using the new positions
	for( i=0; i<NS; i++ ) {
		//Same swimmer
		smonoForce_sameSwimmer( a,SS,(swimmers+i),springType );
		for( d=0; d<DIM; d++ ) (swimmers+i)->H.A[d] += a[d]*SS.iheadM;
		for( d=0; d<DIM; d++ ) (swimmers+i)->M.A[d] -= a[d]*SS.imiddM;

		smonoForce_AlignmentPot(a, CL,SS,(swimmers + i));
		// for( d=0; d<DIM; d++ ) (swimmers+i)->H.A[d] += a[d]*SS.iheadM;
		// for( d=0; d<DIM; d++ ) (swimmers+i)->M.A[d] -= a[d]*SS.imiddM;

		//swimmer alignment potential
		// TODO: replicate lines 461-463 (near exact!) but using 
		//			smonoForce_AlignmentPot. You will also need to pass it the 
		// 			parent MPCD cell as CL

		//All other swimmers
		for( j=i+1; j<NS; j++ ) {
			//Head-head
			smonoForce_differentSwimmers( a,SS,(swimmers+i)->H,(swimmers+j)->H );
			for( d=0; d<DIM; d++ ) (swimmers+i)->H.A[d] += a[d]*SS.iheadM;
			for( d=0; d<DIM; d++ ) (swimmers+j)->H.A[d] -= a[d]*SS.iheadM;
			//Head-middle
			smonoForce_differentSwimmers( a,SS,(swimmers+i)->H,(swimmers+j)->M );
			for( d=0; d<DIM; d++ ) (swimmers+i)->H.A[d] += a[d]*SS.iheadM;
			for( d=0; d<DIM; d++ ) (swimmers+j)->M.A[d] -= a[d]*SS.imiddM;
			//Middle-head
			smonoForce_differentSwimmers( a,SS,(swimmers+i)->M,(swimmers+j)->H );
			for( d=0; d<DIM; d++ ) (swimmers+i)->M.A[d] += a[d]*SS.imiddM;
			for( d=0; d<DIM; d++ ) (swimmers+j)->H.A[d] -= a[d]*SS.iheadM;
			//Middle-middle
			smonoForce_differentSwimmers( a,SS,(swimmers+i)->M,(swimmers+j)->M );
			for( d=0; d<DIM; d++ ) (swimmers+i)->M.A[d] += a[d]*SS.imiddM;
			for( d=0; d<DIM; d++ ) (swimmers+j)->M.A[d] -= a[d]*SS.imiddM;
		}
	}

	//Verlet step 4
	//Update velocity
	for( i=0; i<NS; i++ ) for( d=0; d<DIM; d++ ) {
		(swimmers+i)->H.V[d] += 0.5*dt*(swimmers+i)->H.A[d];
		(swimmers+i)->M.V[d] += 0.5*dt*(swimmers+i)->M.A[d];
	}
}
void smonoDist( double r[],double *dr,smono m1, smono m2 ) {
	//Calculate the distance between two swimmer monomers
	int d;
	for( d=0; d<DIM; d++ ) r[d] = m1.Q[d]-m2.Q[d];
	//Apply any necessary periodic BCs
	swimmerPBC_dr( r );
	*dr = length( r,DIM );
}
double smonoForceMag_sameSwimmer( double dr,specSwimmer SS,swimmer *s,int springType ) {
	//Calculate the forces between two monomers in the same swimmer
	//Use the instantaneous
	double fWCA,fSPRING;
	//WCA
	dr *= (s->isig);		//scale dr for WCA
	fWCA = swimmerWCA( dr,SS.eps );
	//FENE
	dr *= (s->sig)*(s->iro);		//scale dr for FENE
	if(springType==FENESPRING) fSPRING = swimmerFENE( dr,SS.k );
	else if(springType==HOOKESPRING) fSPRING = swimmerHookean( dr,SS.k );
	else fSPRING=0.0;
	return fWCA+fSPRING;
}
void smonoForce_sameSwimmer( double a[],specSwimmer SS,swimmer *s,int springType ) {
	//Calculate the forces between two monomers in the same swimmer
	double dr,f,r[DIM];
	int d;
	//Separation
	smonoDist( r,&dr,s->H,s->M );
	//Forces
	f=smonoForceMag_sameSwimmer( dr,SS,s,springType );
	#ifdef DBG
			if( DBUG == DBGSWIMMER || DBUG == DBGSWIMMERDEETS ) if( fabs(f)*SS.iheadM > 1000.0 || fabs(f)*SS.imiddM > 1000.0) {
				printf("Warning: Swimmer MD forces large\n\tf=%e\n\tSeparation=%lf\n\tRelative position:",f,dr);
				pvec(r,DIM);
			}
	#endif
	//a=acceleration time mass (still a force)
	for( d=0; d<DIM; d++ ) a[d] = f*r[d];
}
void smonoForce_AlignmentPot(double a[], cell ***CL, specSwimmer SS, swimmer *s){
	// //Define the potential type
	// // 1 = Angular Harmonic, 2 = Cosine Harmonic, 3 = Cosine Expansion
	// int x = 0, y = 0, z = 0; //Position in the LC
	// const int pot = 1;
	// const int potConst = 1; //Spring relaxation constant
	// double bacOriLen, forceCoefs, bacOriNorm[DIM], localDir[DIM], theta;

 	// //Calculate norm of r (bacteria body -> head vector)
	
 	// smonoDist(bacOriNorm,&bacOriLen,s->H,s->M);  //Calcualtes distance between body to head
 	// swimmerOri(bacOriNorm,s); //Returns normal vector between body to head

	// x = (int)s->M.Q[0];    //Assigns the bacteria body position as local LC indicies (x,y,z)
	// if (DIM > 1) y = (int)s->M.Q[1]; //Acconts for multiple dimensions 
	// if (DIM >2) z = (int)s->M.Q[2];

	// for (int i = 0; i <DIM; i++) localDir[i] = CL[x][y][z].DIR[i]; // local director at the bacteria's position (a,b,c)


	// // Calculate the angle between the bacteria and the director
 	// theta = acos(dotprod(bacOriNorm, localDir, DIM));

	// //Calculate the force with according to different potentials (pot = 1,2,3)
	// if (pot == 1){ //Angular Harmonic Force 
	// 	forceCoefs = potConst*theta/(sin(theta)*bacOriLen); //Force coefficients independant of componant
	// 	for (int i = 0; i <DIM; i++){
	// 		a[i] = forceCoefs*(cos(theta)*bacOriNorm[i] - localDir[i]);
	// 	}

	// } else if (pot == 2){ //Cosine Harmonic force
	// 	forceCoefs = -2*potConst*cos(theta)/bacOriLen; //Force coefficients 
	// 	for (int i = 0; i <DIM; i++){
	// 		a[i] = forceCoefs*(cos(theta)*bacOriNorm[i] - localDir[i]);
	// 	}

	// } else if (pot == 3){ //Cosine Expansion
	// 	forceCoefs = potConst/bacOriLen; //Force coefficients
	// 	for (int i = 0; i <DIM; i++){
	// 		a[i] = forceCoefs*(cos(theta)*bacOriNorm[i] - localDir[i]);
	// 	}
	// }


// 	/*
// 		Args:
// 			a[] - Acceleration vector to return
// 			*CL - A pointer to the PARENT cell of the swimmer
// 			SS - List of species information on swimmers
// 			*s - A pointer to the particular swimmer instance

// 		Note to Tom:
// 			This is where you want to implement your bending potential. You want
// 			to pass all the things you need to compute it, as a paremeter, down
// 			to this function. 
// 			CL will contain the director of the parent cell.
// 			You will have two hardcoded variables:
// 				- A hardcoded switch on which potential to us
// 				- Bending potential constant
// 			At the end of this function, argument a MUST contain the acceleration
// 			due to your aligment potential


//			For the crystal: magTorque_CL( CL[x][y][z],SP,(double)dt,bacOriNorm* magnitude of the torque ); 

// 	*/
// 	// TODO: write this
}
double smonoForceMag_differentSwimmers( double dr,specSwimmer SS ) {
	//Calculate the forces between two monomers in different swimmers (no FENE)
	//Other swimmers always see the "true" size --- never the shrunk size of tumbling
	//WCA
	dr *= SS.isig;		//scale dr for WCA
	return swimmerWCA( dr,SS.eps );
}
void smonoForce_differentSwimmers( double a[],specSwimmer SS,smono s1,smono s2 ) {
	//Calculate the forces between two monomers in the same swimmer
	double dr,f,r[DIM];
	int d;
	//Separation
	smonoDist( r,&dr,s1,s2 );
	//Forces
	f=smonoForceMag_differentSwimmers( dr,SS );
	#ifdef DBG
			if( DBUG == DBGSWIMMER || DBUG == DBGSWIMMERDEETS ) if( fabs(f)*SS.iheadM > 1000.0 || fabs(f)*SS.imiddM > 1000.0) {
				printf("Warning: Swimmer MD forces large\n\tf=%e\n\tSeparation=%lf\n\tRelative position:",f,dr);
				pvec(r,DIM);
			}
	#endif
	//a=acceleration time mass (still a force)
	for( d=0; d<DIM; d++ ) a[d] = f*r[d];
}

void integrateSwimmers( specSwimmer SS,swimmer swimmers[],bc WALL[],int stepsMD,double timeStep,double MAG[],int springType, cell ***CL) {
/*
    Integrate the motion of the swimmer
*/
	int t,i;
	double dt = timeStep/(double)stepsMD;

	if( SS.TYPE!=DUMBBELL_FIXED ) {
		#ifdef DBG
			if( DBUG == DBGSWIMMER ) {
				printf( "\tPre-integration:\n" );
				for( i=0; i<NS; i++ ) {
					printf( "\t\tS%d:\n",i );
					swcoord(swimmers[i]);
				}
			}
		#endif
		//Dumbbell swimmer
		//Currently swimmers DO NOT have steric interactions BETWEEN swimmers
		for( t=0; t<stepsMD; t++ ) {
			#ifdef DBG
				if( DBUG == DBGSWIMMERDEETS ) {
					printf( "\tMD t=%d:\n",t );
					for( i=0; i<NS; i++ ) {
						printf( "\tS%d:\n",i );
						swcoord(swimmers[i]);
					}
				}
			#endif
			//Apply magnetic field to magnetotaxic swimmer
			#ifdef DBG
				if( DBUG == DBGMAG || DBUG == DBGSWIMMERDEETS ) printf( "\tApply magnetic field.\n" );
			#endif
			for( i=0; i<NS; i++ ) swimmerMagTorque( SS,(swimmers+i),dt,MAG );
			#ifdef DBG
				if( DBUG == DBGSWIMMERDEETS ) for( i=0; i<NS; i++ ) {
					printf( "\tS%d:\n",i );
					swcoord(swimmers[i]);
				}
			#endif
			//Verlet time step
	    #ifdef DBG
	      if( DBUG == DBGSWIMMER || DBUG == DBGSWIMMERDEETS ) printf( "\tVelocity Verlet.\n" );
	    #endif
			// Velocity Verlet integration
			if(SS.TYPE==DUMBBELL_EXVOL) swimmerVerlet_all( SS,swimmers,dt,springType, CL );
			else for( i=0; i<NS; i++ ) swimmerVerlet_nonInteracting( SS,(swimmers+i),dt,springType, CL );
			#ifdef DBG
				if( DBUG == DBGSWIMMERDEETS ) for( i=0; i<NS; i++ ) {
					printf( "\tS%d:\n",i );
					swcoord(swimmers[i]);
				}
			#endif
	    // Swimmer BC
			for( i=0; i<NS; i++ ) {
		    //Head
		    #ifdef DBG
		      if( DBUG == DBGSWIMMER || DBUG == DBGSWIMMERDEETS ) {
		        printf( "\tBCs S%d:\n",i );
		        printf( "\t\tHead bead BCs.\n" );
		      }
		    #endif
		    swimmer_BCcollision( &(swimmers[i].H),WALL,SS,dt );
		    //Middle
		    #ifdef DBG
		      if( DBUG == DBGSWIMMER || DBUG == DBGSWIMMERDEETS ) printf( "\t\tMiddle bead BCs.\n" );
		    #endif
		    swimmer_BCcollision( &(swimmers[i].M),WALL,SS,dt );
				#ifdef DBG
					if( DBUG == DBGSWIMMERDEETS ) for( i=0; i<NS; i++ ) {
						printf( "\tS%d:\n",i );
						swcoord(swimmers[i]);
					}
				#endif
			}
			//If held near wall then set z=1 (in 3D) or y=1 (in 2D)
			if( SS.TYPE==DUMBBELL_NEARWALL ) for( i=0; i<NS; i++ ) {
				(swimmers+i)->H.Q[DIM-1]=SS.fixDist;
				(swimmers+i)->H.V[DIM-1]=0.0;
				(swimmers+i)->M.Q[DIM-1]=SS.fixDist;
				(swimmers+i)->M.V[DIM-1]=0.0;
			}
		}
		#ifdef DBG
			if( DBUG == DBGSWIMMER ) {
				printf( "\tPost-integration:\n" );
				for( i=0; i<NS; i++ ) {
					printf( "\t\tS%d:\n",i );
					swcoord(swimmers[i]);
				}
			}
		#endif
	}
}
void swimmerMagTorque( specSwimmer SS,swimmer *sw,double dt,double MAG[] ) {
	/*
	    Apply the magnetic torque to all the swimmers
	*/
	int d;
	double magT,R;
	double force[DIM],r[DIM],u[_3D],n[_3D],mT[_3D];

	for( d=0; d<_3D; d++) u[d]=0.0;

	//Vector from middle to head
	for( d=0; d<DIM; d++ ) r[d] = sw->H.Q[d]-sw->M.Q[d];
	swimmerPBC_dr( r );
	//Normalize to get the direction of the magnetic moment of the magnetotaxic swimmer
	normCopy( r,u,DIM );
	//Calculate the cross product of u (proportional to the magnetic moment) and magnetic field (which is then in the direction of the torque)
	crossprod( u,MAG,mT );
	//Find the torque by multiplying by the strength of the magnetic moment
	for( d=0; d<_3D; d++ ) mT[d] *= SS.MAGMOM;
	#ifdef DBG
		if( DBUG == DBGMAG || DBUG == DBGSWIMMERDEETS ) {
			printf( "Swimmer orientation: " );
			pvec( u,_3D );
			printf( "Magnetic field: ");
			pvec( MAG,_3D );
			printf( "Torque: " );
			pvec( mT,_3D );
		}
	#endif
	//Calculate the direction of the force on the head
	crossprod( mT,u,n );
	magT=length( n,_3D );
	for( d=0; d<_3D; d++) n[d]/=magT;
	//Calculate the distance from the CM to the head (assuming same mass)
	R=0.5*length( r,DIM );
	//Magnitude of the torque
	magT=length( mT,_3D );
	#ifdef DBG
		if( DBUG == DBGMAG || DBUG == DBGSWIMMERDEETS ) {
			printf( "Mag torque: %lf\n",magT );
		}
	#endif
	//Calculate the force on the head and the body needed to produce the torque
	if(magT>0.0) for( d=0; d<DIM; d++ ) force[d] = 0.5*magT*n[d]/R;
	else for( d=0; d<DIM; d++ ) force[d] = 0.0;
	#ifdef DBG
		if( DBUG == DBGMAG || DBUG == DBGSWIMMERDEETS  ) {
			printf( "Force: " );
			pvec( force,DIM );
			printf( "Percent Diff Head:" );
			if(DIM==_3D) printf( " (%lf, %lf, %lf)\n", 100*force[0]*SS.iheadM*dt/sw->H.V[0], 100*force[1]*SS.iheadM*dt/sw->H.V[1], 100*force[2]*SS.iheadM*dt/sw->H.V[2] );
			else if(DIM==_2D) printf( " (%lf, %lf)\n", 100*force[0]*SS.iheadM*dt/sw->H.V[0], 100*force[1]*SS.iheadM*dt/sw->H.V[1] );
			printf( "Percent Diff Midd:" );
			if(DIM==_3D) printf( " (%lf, %lf, %lf)\n", 100*force[0]*SS.imiddM*dt/sw->M.V[0], 100*force[1]*SS.imiddM*dt/sw->M.V[1], 100*force[2]*SS.imiddM*dt/sw->M.V[2] );
			else if(DIM==_2D) printf( " (%lf, %lf)\n", 100*force[0]*SS.imiddM*dt/sw->M.V[0], 100*force[1]*SS.imiddM*dt/sw->M.V[1] );
		}
	#endif
	//Apply force/impulse to swimmer's head
	for( d=0; d<DIM; d++ ) sw->H.V[d] += force[d]*SS.iheadM*dt;
	//Apply equal but opposite force/impulse to swimmer's body
	for( d=0; d<DIM; d++ ) sw->M.V[d] -= force[d]*SS.imiddM*dt;
}
void allSwimmersMagTorque( specSwimmer SS,swimmer swimmers[],double timeStep,int stepsMD,double MAG[] ) {
	/*
	    Apply the magnetic torque to all the swimmers
	*/
	int i,t,d;
	double magT,R,dt;
	double force[DIM],r[DIM],u[_3D],n[_3D],mT[_3D];

	for( d=0; d<_3D; d++) u[d]=0.0;
	dt=timeStep/(double)stepsMD;

	for( i=0; i<NS; i++ ) for( t=0; t<stepsMD; t++ ) {
		//Vector from middle to head
		for( d=0; d<DIM; d++ ) r[d] = (swimmers+i)->H.Q[d]-(swimmers+i)->M.Q[d];
		swimmerPBC_dr( r );
		//Normalize to get the direction of the magnetic moment of the magnetotaxic swimmer
		normCopy( r,u,DIM );
		//Calculate the cross product of u and magnetic field (which is in the direction of the torque)
		crossprod( u,MAG,mT );
		//Find the torque by multiplying by the strength of the magnetic moment
		for( d=0; d<_3D; d++ ) mT[d] *= SS.MAGMOM;
		#ifdef DBG
			if( DBUG == DBGMAG ) {
				printf( "Swimmer %d orientation sub-timestep %d: ",i,t );
				pvec( u,_3D );
				printf( "Magnetic field: ");
				pvec( MAG,_3D );
				printf( "Torque: " );
				pvec( mT,_3D );
			}
		#endif
		//Calculate the direction of the force on the head
		crossprod( mT,u,n );
		magT=length( n,_3D );
		for( d=0; d<_3D; d++) n[d]/=magT;
		//Calculate the distance from the CM to the head
		R=0.5*length( r,DIM );
		//Magnitude of the torque
		magT=length( mT,_3D );
		//Calculate the force on the head and the body needed to produce the torque
		if(magT>0.0) for( d=0; d<DIM; d++ ) force[d] = 0.5*magT*n[d]/R;
		else for( d=0; d<DIM; d++ ) force[d] = 0.0;
		#ifdef DBG
			if( DBUG == DBGMAG ) {
				printf( "Force: " );
				pvec( force,DIM );
			}
		#endif
		//Apply force/impulse to swimmer's head
		for( d=0; d<DIM; d++ ) (swimmers+i)->H.V[d] += force[d]*SS.iheadM*dt;
		//Apply equal but opposite force/impulse to swimmer's body
		for( d=0; d<DIM; d++ ) (swimmers+i)->M.V[d] -= force[d]*SS.imiddM*dt;
	}
}
void swimmerDipole( specSwimmer SS,swimmer swimmers[],cell ***CL,spec SP[],double timeStep,particleMPC *SRDparticles,bc WALL[],simptr simMD ) {
/*
    Apply both the force dipole and the rotlet-torque dipole to each swimmer
*/
	int i,d;
	int aH,bH,cH,aT,bT,cT;
	double r[DIM],QT[_3D];	// Position of tail
	double SHIFT[_3D];			// Vector to positively shift all components of the simulation

	for( d=0; d<_3D; d++ ) QT[d] = 0.0;
	for( d=0; d<_3D; d++ ) SHIFT[d] = 0.0;

	for( i=0; i<NS; i++ ) if( (swimmers+i)->RT==RUNNING ) {

		/* ****************************************** */
		/* ************** FORCE DIPOLE ************** */
		/* ****************************************** */
		#ifdef DBG
			if( DBUG == DBGSWIMMERDEETS ) {
				printf( "Apply Force Dipole to Swimmer %d\n",i);
			}
		#endif
		if( !feq(fabs(SS.FS),0.0) )	swimmerForceDipole( SS,(swimmers+i),CL,SP,timeStep );

		/* ****************************************** */
		/* ********** ROTLET/TORQUE DIPOLE ********** */
		/* ****************************************** */
		if( DIM>=_3D && !feq(fabs(SS.TS),0.0) )	{

			#ifdef DBG
				if( DBUG == DBGSWIMMERDEETS || DBUG == DBGSWIMMERTORQUE ) {
					printf( "Apply Torque Dipole to Swimmer %d\n",i);
				}
			#endif
			/* ****************************************** */
			/* ********** GRID SHIFT & BINNING ********** */
			/* ****************************************** */
			#ifdef DBG
				if( DBUG == DBGSWIMMERTORQUE ) printf( "\tRotlet Shift Grid.\n" );
			#endif
			//Generate shift
			//Find which cell the head is in
			aH = (int)(swimmers+i)->H.Q[0];
			bH = (int)(swimmers+i)->H.Q[1];
			cH = (int)(swimmers+i)->H.Q[2];
			SHIFT[0] = aH+0.5-(swimmers+i)->H.Q[0];
			if( DIM >= _2D ) SHIFT[1] = bH+0.5-(swimmers+i)->H.Q[1];
			if( DIM >= _3D ) SHIFT[2] = cH+0.5-(swimmers+i)->H.Q[2];
			#ifdef DBG
				if( DBUG == DBGSWIMMERTORQUE ) {
					printf( "\tHead =" );
					pvec((swimmers+i)->H.Q,DIM);
					printf( "\tShift =" );
					pvec(SHIFT,DIM);
					printf( "\tRotlet Shift Grid.\n" );
				}
			#endif
			//Shift entire system by SHIFT
			gridShift_all( SHIFT,0,SRDparticles,WALL,simMD,swimmers,noMD );
			// Bin
			#ifdef DBG
				if( DBUG == DBGSWIMMERTORQUE ) printf( "\tRotlet Bin Particles.\n" );
			#endif
			// Bin SRD particles
			bin( CL,SP,WALL,1.0,ISOF,1 );
			// Bin swimmer monomers
			binSwimmers( CL,1 );
			#ifdef DBG
				if( DBUG == DBGSWIMMERTORQUE ) {
					printf( "\tShifted Head =" );
					pvec((swimmers+i)->H.Q,DIM);
				}
			#endif
			// Find which cell the tail is now in
			for( d=0; d<DIM; d++ ) r[d] = (swimmers+i)->H.Q[d] - (swimmers+i)->M.Q[d];
			swimmerPBC_dr( r );
			for( d=0; d<DIM; d++ ) QT[d] = (swimmers+i)->M.Q[d] - SS.DS*r[d];
			//Apply PBCs to this tail position (just do PBCs)
			for( d=0; d<DIM; d++ ) {
				if( QT[d] < 0.0 ) QT[d] += XYZ[d];
				else if( QT[d] >= XYZ[d] ) QT[d] -= XYZ[d];
			}
			aT = (int)QT[0];
			bT = (int)QT[1];
			cT = (int)QT[2];
			// Calculate new cell properties for head/tail cells
			CL[aH][bH][cH].POP = localPOP( CL[aH][bH][cH] );
			CL[aH][bH][cH].MASS = localMASS( CL[aH][bH][cH],SP,SS );
			localCM( &CL[aH][bH][cH],SP,SS );
			CL[aT][bT][cT].POP = localPOP( CL[aT][bT][cT] );
			CL[aT][bT][cT].MASS = localMASS( CL[aT][bT][cT],SP,SS );
			localCM( &CL[aT][bT][cT],SP,SS );
			#ifdef DBG
				if( DBUG == DBGSWIMMERTORQUE ) {
					printf( "\tNumber of SRD particles in shifted head pos = %d\n",CL[aH][bH][cH].POP );
				}
			#endif





			/* ****************************************** */
			/* ************* ROTLET DIPOLE ************** */
			/* ****************************************** */
			// Apply the rotlet to a cell centred on the swimmers head
			swimmerRotletDipole( SS,(swimmers+i),CL,SP,timeStep );



			/* ****************************************** */
			/* ************ GRID SHIFT BACK ************* */
			/* ****************************************** */
			#ifdef DBG
				if( DBUG == DBGSWIMMERTORQUE ) printf( "\tRotlet Shift Grid Back.\n" );
			#endif
			//Shift entire system by back
			gridShift_all( SHIFT,1,SRDparticles,WALL,simMD,swimmers,noMD );
			/* ****************************************** */
			/* ***************** RE-BIN ***************** */
			/* ****************************************** */
			#ifdef DBG
				if( DBUG >= DBGTITLE ) printf( "\tRe-bin Particles.\n" );
			#endif
			// Bin SRD particles
			bin( CL,SP,WALL,1.0,ISOF,0 );
			// Bin swimmer monomers
			binSwimmers( CL,0 );
			// Re-calculate old cell properties for head/tail cells
			CL[aH][bH][cH].POP = localPOP( CL[aH][bH][cH] );
			CL[aH][bH][cH].MASS = localMASS( CL[aH][bH][cH],SP,SS );
			localCM( &CL[aH][bH][cH],SP,SS );
			CL[aT][bT][cT].POP = localPOP( CL[aT][bT][cT] );
			CL[aT][bT][cT].MASS = localMASS( CL[aT][bT][cT],SP,SS );
			localCM( &CL[aT][bT][cT],SP,SS );
		}
	}
}
void swimmerForceDipole( specSwimmer SS,swimmer *sw,cell ***CL,spec SP[],double timeStep ) {
/*
    Set the swimming speed and the propulsion force on the fluid
    This is where the "invisible" tail comes in.
*/
  int a=0,b=0,c=0,d=0;
	double r[DIM],n[DIM],acc[DIM],QT[_3D];
	double m;
	particleMPC *pMPC;		//Temporary pointer to MPC particles
	particleMD *pMD;	//Temporary pointer to MD particles
	smono *pSW;			//Temporary pointer to swimmer monomers

  for( d=0; d<_3D; d++ ) QT[d] = 0.0;
  for( d=0; d<DIM; d++ ) {
    r[d] = 0.0;
    n[d] = 0.0;
    acc[d] = 0.0;
  }

  //Vector from middle to head
	for( d=0; d<DIM; d++ ) r[d] = sw->H.Q[d] - sw->M.Q[d];
	swimmerPBC_dr( r );
	normCopy( r,n,DIM );
	//Apply force to swimmer's head (half of force on head/half on middle)
	for( d=0; d<DIM; d++ ) acc[d] = 0.5*SS.FS*n[d]*SS.iheadM;
	#ifdef DBG
		if( DBUG == DBGSWIMMERDEETS ) {
			printf( "\t\tHead mass=%d and acc = ", SS.headM);
			pvec( acc,DIM );
		}
	#endif
	for( d=0; d<DIM; d++ ) sw->H.V[d] += acc[d]*timeStep;
	//Apply force to swimmer's middle body (half of force on head/half on middle)
	for( d=0; d<DIM; d++ ) acc[d] = 0.5*SS.FS*n[d]*SS.imiddM;
	#ifdef DBG
		if( DBUG == DBGSWIMMERDEETS ) {
			printf( "\t\tBody mass=%d and acc = ", SS.headM);
			pvec( acc,DIM );
		}
	#endif
	for( d=0; d<DIM; d++ ) sw->M.V[d] += acc[d]*timeStep;
  //The force on the "invisible tail" balances the force needed to move the head propulsion
	if(SS.TYPE!=DUMBBELL_MONOF) {
		//Find the position of the "phantom" tail around the CM position
		for( d=0; d<DIM; d++ ) QT[d] = sw->M.Q[d] - SS.DS*r[d];
		//Apply PBCs to this tail position (just do PBCs)
		for( d=0; d<DIM; d++ ) {
			if( QT[d] < 0.0 ) QT[d] += XYZ[d];
			else if( QT[d] >= XYZ[d] ) QT[d] -= XYZ[d];
		}
		//Find which cell the tail is in
		a = (int)QT[0];
		b = (int)QT[1];
		c = (int)QT[2];
		//Apply the force to the fluid in the tail's cell
    if( CL[a][b][c].MASS>0.0 ) {
			//Calculate the acceleration per particle in the cell (balances force on head and middle body)
			for( d=0; d<DIM; d++ ) acc[d] = -SS.FS*n[d] / CL[a][b][c].MASS;
			#ifdef DBG
				if( DBUG == DBGSWIMMERDEETS ) {
					printf( "\t\tCell [%d,%d,%d] mass=%lf and acc = ", a,b,c,CL[a][b][c].MASS);
					pvec( acc,DIM );
				}
			#endif
  		// Apply the force to the SRD particles
  		if( CL[a][b][c].pp != NULL ) {
  			pMPC = CL[a][b][c].pp;
  			while(pMPC != NULL) {
					m=SP[pMPC->SPID].MASS;
  				for( d=0; d<DIM; d++ ) pMPC->V[d] += acc[d]*timeStep*m;
  				//Increment link in list
  				pMPC = pMPC->next;
  			}
			}
			// MD particles
			if(MDmode==MDinMPC) if( CL[a][b][c].MDpp!=NULL) {
				pMD = CL[a][b][c].MDpp;
				while( pMD!=NULL ) {
					m = pMD->mass;
					pMD->vx += acc[0]*timeStep*m;
					if( DIM >= _2D ) pMD->vy += acc[1]*timeStep*m;
					if( DIM >= _3D ) pMD->vz += acc[2]*timeStep*m;
					//Increment link in list
					pMD = pMD->nextSRD;
				}
			}
			// Swimmer particles
			if( CL[a][b][c].sp!=NULL) {
				pSW = CL[a][b][c].sp;
				while( pSW!=NULL ) {
					if( pSW->HorM ) m = (double) SS.middM;
					else m = (double) SS.headM;
					for( d=0; d<DIM; d++ ) pSW->V[d] += acc[d]*timeStep*m;
					//Increment link in list
					pSW = pSW->next;
				}
			}
		}
	}
}
void swimmerRotletDipole( specSwimmer SS,swimmer *sw,cell ***CL,spec SP[],double timeStep ) {
/*
    Set the swimmer's rotations and the torque on the fluid
    This is also where the "invisible" tail comes in.
*/
	int a=0,b=0,c=0,d=0;
	double q_sw[_3D];													//Position of the swimmers' head or tail
	double r_mh[DIM],n_mh[DIM]; 							//Vectors betweem middle/body and head
	double r_cm[DIM],momI; 										//Centre of mass position of fluid in cell and moment of inertia
	double dw;																//Mag of angular velocity on each MPCD particle

	for( d=0; d<_3D; d++ ) q_sw[d] = 0.0;
	for( d=0; d<DIM; d++ ) {
		r_mh[d] = 0.0;
		n_mh[d] = 0.0;
		r_cm[d] = 0.0;
	}
	momI = 0.0;

	//Vector from middle to head
	for( d=0; d<DIM; d++ ) r_mh[d] = sw->H.Q[d] - sw->M.Q[d];
	swimmerPBC_dr( r_mh );
	normCopy( r_mh,n_mh,DIM );
	if(SS.TYPE!=DUMBBELL_MONOF) {

		//Apply the rotlet to the HEAD
		//Save position of head
		for( d=0; d<DIM; d++ ) q_sw[d] = sw->H.Q[d];
		//Apply torque to MPCD particles in cell that includes swimmer's head
		//Find which cell the head is in
		a = (int)q_sw[0];
		b = (int)q_sw[1];
		c = (int)q_sw[2];
		#ifdef DBG
			if( DBUG == DBGSWIMMERTORQUE ) {
				printf( "\tHead cell [%d,%d,%d] mass %lf; head mass %d\n",a,b,c,CL[a][b][c].MASS,SS.headM );
				printf( "\t\tSwimmer direction" );
				pvec(n_mh,DIM);
			}
		#endif
		//Apply the rotational forces to the fluid in the head's cell (need to subtract head's mass)
		if( CL[a][b][c].MASS-SS.headM>0.0 ) {
			// Calculate the centre of mass
			localCM_SRD( CL[a][b][c],SP,r_cm );
			// Calculate the moment of inertia about the centre of mass and direction of swimmer
			momI=localMomInertia_SRD( CL[a][b][c],SP,r_cm,n_mh );
			//Calculate the change in angular speed
			dw = timeStep*SS.TS/momI;
			#ifdef DBG
				if( DBUG == DBGSWIMMERTORQUE ) {
					printf( "\t\tI=%lf, dw=%lf, CM=",momI,dw );
					pvec(r_cm,DIM);
				}
			#endif
			// Apply the torque to the SRD particles
			rotate_CL( CL[a][b][c],SP,r_cm,n_mh,dw );
		}

		//Apply the equal-but-opposite rotlet to the TAIL
		for( d=0; d<DIM; d++ ) r_cm[d] = 0.0;
		//Save position of tail
		for( d=0; d<DIM; d++ ) q_sw[d] = sw->M.Q[d] - SS.DS*r_mh[d];
		//Apply PBCs to this tail position (just do PBCs)
		for( d=0; d<DIM; d++ ) {
			if( q_sw[d] < 0.0 ) q_sw[d] += XYZ[d];
			else if( q_sw[d] >= XYZ[d] ) q_sw[d] -= XYZ[d];
		}
		//Apply torque to MPCD particles in cell that includes swimmer's "invisible tail"
		//Find the position of the "phantom" tail around the CM position
		//Find which cell the tail is in
		a = (int)q_sw[0];
		b = (int)q_sw[1];
		c = (int)q_sw[2];
		#ifdef DBG
			if( DBUG == DBGSWIMMERTORQUE ) {
				printf( "\tTail cell [%d,%d,%d] mass %lf; head mass %d\n",a,b,c,CL[a][b][c].MASS,SS.headM );
			}
		#endif
		//Apply the rotational forces to the fluid in the head's cell
		if( CL[a][b][c].MASS>0.0 ) {
			// Calculate the centre of mass
			localCM_SRD( CL[a][b][c],SP,r_cm );
			// Calculate the moment of inertia about the centre of mass and direction of swimmer
			momI=localMomInertia_SRD( CL[a][b][c],SP,r_cm,n_mh );
			//Calculate the change in angular speed
			dw = -timeStep*SS.TS/momI;
			#ifdef DBG
				if( DBUG == DBGSWIMMERTORQUE ) {
					printf( "\t\tI=%lf, dw=%lf, CM=",momI,dw);
					pvec(r_cm,DIM);
				}
			#endif
			// Apply the torque to the SRD particles
			rotate_CL( CL[a][b][c],SP,r_cm,n_mh,dw );
		}
	}
}
void runTumbleDynamics( specSwimmer *SS,swimmer swimmers[],bc WALL[],int stepsMD,double MAG[],double dt,int RTOUT,FILE *fruntumble, cell ***CL) {
/*
    Stochastically run and tumble.
    This routine checks the number of times since last run/tumble switching event.
    If an event occurs a new random run/tumble time is generated
    If the swimmer tumbles then it's size can shrink (or techniqually grow but this shouldn't occur)
    If it runs then it'shrinkMDSteps size is returned to normal
		Uses a hookean spring during shrinking so that FENE spring isn't broken
*/
  int i,j;
  double dr,ds,dk;			//FENE separation and LJ sigma step size for changing
  int shrinkMDSteps=10*stepsMD;
  double dt_int;
	double n[DIM];				//Swimmer orientation

	#ifdef DBG
		if( DBUG == DBGRUNTUMBLE ) for( i=0; i<NS; i++ ) {
			printf( "\tSwimmer %d. Run/tumble/shrink/extend phase %d. Period=%d. Time=%d\n", i, (swimmers+i)->RT, (swimmers+i)->timeRND, (swimmers+i)->timeCNT );
		}
	#endif
	//Since the routine integrates during shrinking, must breakup integration time
	if( fneq(SS->sizeShrink,1.0) && fneq(SS->springShrink,1.0) ) dt_int=dt/(double)shrinkMDSteps/3.0;
	else if( fneq(SS->sizeShrink,1.0) ) dt_int=dt/(double)shrinkMDSteps/2.0;
	else dt_int=dt/(double)shrinkMDSteps;

	for( i=0; i<NS; i++ ) {
		/* ****************************************** */
		/* *************** Run/Tumble *************** */
		/* ****************************************** */
    //Running phase
    if( (swimmers+i)->RT==RUNNING ){
      //Haven't reached tumbling time yet
      if( (swimmers+i)->timeCNT < (swimmers+i)->timeRND ) (swimmers+i)->timeCNT++;
      //Reached tumbling time --- switch to shrinking phase
      else {
				if( RTOUT>=OUT ) {
					runtumbleout( fruntumble,(swimmers+i) );
					swimmerOri( n,(swimmers+i) );
					for( j=0; j<DIM; j++ ) (swimmers+i)->n0[j]=n[j]; //Reset old orientation to current orientation
				}
        (swimmers+i)->RT = SHRINKING;      //Set in shrinking phase
        (swimmers+i)->timeCNT = 0;        //Reset counter
        #ifdef DBG
          if( DBUG == DBGRUNTUMBLE ) {
            printf( "Run->Shrink mode\n" );
          }
        #endif
      }
    }
    //Tumbling phase
    else if( (swimmers+i)->RT==TUMBLING ){
      //Haven't reached running time yet
      if( (swimmers+i)->timeCNT < (swimmers+i)->timeRND ) (swimmers+i)->timeCNT++;
      //Reached running time --- switch to extending phase
      else {
				if( RTOUT>=OUT ) {
					runtumbleout( fruntumble,(swimmers+i) );
					swimmerOri( n,(swimmers+i) );
					for( j=0; j<DIM; j++ ) (swimmers+i)->n0[j]=n[j]; //Reset old orientation to current orientation
				}
        (swimmers+i)->RT = EXTENDING;      //Set in extending phase
        (swimmers+i)->timeCNT = 0;        //Reset counter
				#ifdef DBG
          if( DBUG == DBGRUNTUMBLE ) {
            printf( "Tumbling->Extend mode\n" );
          }
        #endif
			}
    }
		/* ****************************************** */
		/* ************* Shrink/Extend ************** */
		/* ****************************************** */
		// Shrinking phase
		if( (swimmers+i)->RT==SHRINKING ) {
			//Shrink
			//Shrink from run conformation to tumble conformation
			//Shrink the swimmers slowly (in shrinkMDSteps steps) times more than normal MD
			#ifdef DBG
				if( DBUG == DBGRUNTUMBLE ) {
					printf( "\tInitial ro=%lf, sig=%lf, k=%lf\n",(swimmers+i)->ro,(swimmers+i)->sig,(swimmers+i)->k );
					printf( "\tShrink size coef=%lf\n",SS->sizeShrink );
					printf( "\tShrink spring coef=%lf\n",SS->springShrink );
					printf("\tShrink swimmer into tumbling conformation in %d steps\n",shrinkMDSteps);
				}
			#endif
			if( fneq(SS->sizeShrink,1.0) ) {
				//Shrink size first then shrink spring (factor of 2 because do it in this two step way)
				ds=( SS->sig - (SS->sig)*(SS->sizeShrink) ) / ((double)(shrinkMDSteps * (SS->shrinkTime)));
				for( j=0; j<shrinkMDSteps; j++ ) {
					(swimmers+i)->sig -= ds;
					(swimmers+i)->isig=1.0/((swimmers+i)->sig);
					integrateSwimmers( *SS,swimmers,WALL,shrinkMDSteps,dt_int,MAG,HOOKESPRING,CL );
				}
				//Shrink spring
				dr=( SS->ro - (SS->ro)*(SS->sizeShrink)) / ((double)(shrinkMDSteps * (SS->shrinkTime)));
				for( j=0; j<shrinkMDSteps; j++ ) {
					(swimmers+i)->ro-=dr;
					(swimmers+i)->iro=1.0/((swimmers+i)->ro);
					integrateSwimmers( *SS,swimmers,WALL,shrinkMDSteps,dt_int,MAG,HOOKESPRING, CL );
				}
			}
			if( fneq(SS->springShrink,1.0) ) {
				dk=( SS->k - (SS->k)*(SS->springShrink) ) / ((double)(shrinkMDSteps * (SS->shrinkTime)));
				for( j=0; j<shrinkMDSteps; j++ ) {
					(swimmers+i)->k-=dk;
					integrateSwimmers( *SS,swimmers,WALL,shrinkMDSteps,dt_int,MAG,HOOKESPRING, CL);
				}
			}
			#ifdef DBG
				if( DBUG == DBGRUNTUMBLE ) {
					printf( "\tFinal ro=%lf, sig=%lf, k=%lf\n",(swimmers+i)->ro,(swimmers+i)->sig,(swimmers+i)->k );
				}
			#endif
			//Check if done shrinking
			(swimmers+i)->timeCNT++;        //Iterate counter
			if( (swimmers+i)->timeCNT >= SS->shrinkTime ) {
				(swimmers+i)->RT = TUMBLING;      //Set in tumbling phase
				(swimmers+i)->timeCNT = 0;        //Reset counter
				(swimmers+i)->timeRND = (int) genrand_exp( SS->tumbleTime );    //Generate random integer to tumble for
				#ifdef DBG
					if( DBUG == DBGRUNTUMBLE ) {
						printf( "Shrinking->Tumbling phase\n\tTumbling time=%d\n",(swimmers+i)->timeRND );
					}
				#endif
			}
		}
		//Extending phase
		else if( (swimmers+i)->RT==EXTENDING ) {
				#ifdef DBG
					if( DBUG == DBGRUNTUMBLE ) {
						printf( "Tumble->run\n\tRunning time=%d\n",(swimmers+i)->timeRND );
					}
				#endif
				//Grow from tumble conformation to run conformation
				//Grow the swimmers slowly (in 2*shrinkMDSteps steps) times more than normal MD
				#ifdef DBG
					if( DBUG == DBGRUNTUMBLE ) {
						printf( "\tInitial ro=%lf, sig=%lf, k=%lf\n",(swimmers+i)->ro,(swimmers+i)->sig,(swimmers+i)->k );
						printf( "\tShrink size coef=%lf\n",SS->sizeShrink );
						printf( "\tShrink spring coef=%lf\n",SS->springShrink );
						printf("\tGrow swimmer into running conformation in %d steps\n",shrinkMDSteps);
					}
				#endif
				if( fneq(SS->sizeShrink,1.0) ) {
					dr=( SS->ro-(SS->ro)*(SS->sizeShrink) ) / ((double)(shrinkMDSteps * (SS->shrinkTime)));
					for( j=0; j<shrinkMDSteps; j++ ) {
						(swimmers+i)->ro+=dr;
						(swimmers+i)->iro=1.0/((swimmers+i)->ro);
						integrateSwimmers( *SS,swimmers,WALL,shrinkMDSteps,dt_int,MAG,HOOKESPRING, CL );
					}
					ds=( SS->sig-(SS->sig)*(SS->sizeShrink) ) / ((double)(shrinkMDSteps * (SS->shrinkTime)));
					for( j=0; j<shrinkMDSteps; j++ ) {
						(swimmers+i)->sig+=ds;
						(swimmers+i)->isig=1.0/((swimmers+i)->sig);
						integrateSwimmers( *SS,swimmers,WALL,shrinkMDSteps,dt_int,MAG,HOOKESPRING,CL );
					}
				}
				if( fneq(SS->springShrink,1.0) ) {
					dk=( SS->k - (SS->k)*(SS->springShrink) ) / ((double)(shrinkMDSteps * (SS->shrinkTime)));
					for( j=0; j<shrinkMDSteps; j++ ) {
						(swimmers+i)->k+=dk;
						integrateSwimmers( *SS,swimmers,WALL,shrinkMDSteps,dt_int,MAG,HOOKESPRING,CL);
					}
				}
				#ifdef DBG
					if( DBUG == DBGRUNTUMBLE ) {
						printf( "\tFinal ro=%lf, sig=%lf, k=%lf\n",(swimmers+i)->ro,(swimmers+i)->sig,(swimmers+i)->k );
					}
				#endif
				//Check if done extending
				(swimmers+i)->timeCNT++;        //Iterate counter
				if( (swimmers+i)->timeCNT >= SS->shrinkTime ) {
					(swimmers+i)->RT = RUNNING;      //Set in running phase
					(swimmers+i)->timeCNT = 0;        //Reset counter
					(swimmers+i)->timeRND = (int) genrand_exp( SS->runTime );    //Generate random integer to run for
					#ifdef DBG
						if( DBUG == DBGRUNTUMBLE ) {
							printf( "Extending->Running phase\n\tRunning time=%d\n",(swimmers+i)->timeRND );
						}
					#endif
				}
		}
  }
}
void bininSwimmers( specSwimmer SS,swimmer swimmers[],cell ***CL ) {
/*
   This function does the initial binning of the
   swimmers after they have been first initialized.
   It is different from bin in that it uses the
   actual array of particleMPCs rather than the array
   of pointers to particleMPCs.
*/
	int i,a,b,c;
	//Bin Particles
	for( i=0; i<NS; i++ ){
		//Truncate the coordinate to see what cell the head falls in
		a = (int)swimmers[i].H.Q[0];
		b = (int)swimmers[i].H.Q[1];
		c = (int)swimmers[i].H.Q[2];
		addlinkSwimmer( &CL[a][b][c],&(swimmers[i].H) );
		//Truncate the coordinate to see what cell the middle falls in
		a = (int)swimmers[i].M.Q[0];
		b = (int)swimmers[i].M.Q[1];
		c = (int)swimmers[i].M.Q[2];
		addlinkSwimmer( &CL[a][b][c],&(swimmers[i].M) );
	}
}
void binSwimmers( cell ***CL,int shifted ) {
/*
   This function bins the swimmers i.e. it places
   a pointer to the particleMPC in the appropriate new
   list and removes it from it's old list
*/
	int i,j,k,a,b,c;
	smono *cp;	//Pointer to current item in list
	smono *tp;	//Temporary pointer
	//Search each cell for swimmers that have left the cell
	for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
		cp = CL[i][j][k].sp;
		while(cp != NULL) {
			//Truncate the coordinate to see what cell the particleMPC falls in
			a = (int)cp->Q[0];
			b = (int)cp->Q[1];
			c = (int)cp->Q[2];
			//Save the next address
			tp = cp->next;
			//Apply the periodic BC in binning if the axes have been flagged
			if( shifted ) {
				if( XYZPBC[0] && a==XYZ[0]) a=0;
				if( XYZPBC[1] && b==XYZ[1]) b=0;
				if( XYZPBC[2] && c==XYZ[2]) c=0;
			}

			//Make sure particleMPC didn't escape
			if( a>XYZ[0] || b>XYZ[1] || c>XYZ[2] || a<0 || b<0 || c<0 ){
				printf( "Error:\tSwimmer monomer escaped from cell [%i,%i,%i].\n",a,b,c );
 				exit( 0 );
			}
			//Remove from old list and add to new
			else if( a != i || b != j || c != k ){
				removelinkSwimmer( cp,&CL[i][j][k] );
				addlinkSwimmer( &CL[a][b][c],cp );
			}
			//Return attention to the list under consideration
			cp = tp;
		}
	}
}
void addlinkSwimmer( cell *CL,smono *s ) {
/*
	This routine adds a link to the end of the list
*/
	smono *tp;	//Temporary pointer to swimmer monomer

	if( CL->sp == NULL ) {
		CL->sp = s;
		s->previous = NULL;
	}
	else {
		tp = CL->sp;
		//Find the end of the list
		while( tp->next != NULL ) tp = tp->next;
		//Once the end is found, add the smono to the list
		tp->next = s;
		//Point the smono back at the last link
		s->previous = tp;
	}
	//Swimmer monomer is now at the end of the list so set next to null
	s->next = NULL;
}
void removelinkSwimmer( smono *current,cell *CL ) {
/*
	This routine removes a link from
	a list and relinks the list
*/
	//Point the next link back at the previous link (unless last link)
	if( current->next != NULL ) current->next->previous = current->previous;
	//Point the previous link at the next link (unless first link)
	if( current->previous != NULL ) current->previous->next = current->next;
	//If the current link is the 1st link link the cell to the new first
	else {
		CL->sp = current->next;
	}
	//Remove the current (this is kinda redundant, no?)
	current->previous = NULL;
	current->next = NULL;
}
