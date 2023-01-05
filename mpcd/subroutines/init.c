# include<stdio.h>
# include<math.h>
# include<time.h>
# include<string.h>
# include<stdlib.h>

# include "../headers/definitions.h"
# include "../headers/SRDclss.h"
# include "../headers/globals.h"
# include "../headers/rand.h"
# include "../headers/pout.h"
# include "../headers/mtools.h"
# include "../headers/ctools.h"
# include "../headers/pout.h"
# include "../headers/bc.h"
# include "../headers/mpc.h"
# include "../headers/therm.h"
# include "../headers/lc.h"
# include "../headers/swimmers.h"

# include "../../md/mdtypes.h"
# include "../../md/mdsrd.h"
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* *************** OPEN FILES *************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
void openBasic( FILE **fout,char dir[],char filestring[],char fileextension[] ) {
	char filename[200];

	strcpy( filename,dir );
	strcat( filename,filestring );
	strcat( filename,fileextension );
	*fout = fopen(filename,"w+" );
	if( !*fout ){ // file couldn't be opened
		printf( "Error: File '%s'could not be opened.\n",filename );
		exit( 1 );
	}
	outheader( *fout,0 );
}
void openCheckpoint( FILE **fout,char dir[] ) {
	char filename[200];
	char filechckpoint[]="checkpoint";
	char fileextension[]=".dat";

	strcpy( filename,dir );
	strcat( filename,filechckpoint );
	strcat( filename,fileextension );
	*fout = fopen(filename,"w" );
	if( !*fout ){ // file couldn't be opened
		printf( "Error: File '%s'could not be opened.\n",filename );
		exit( 1 );
	}
}
void opendetails( int i,FILE *fdetail[],char dir[],char fileprefix[],char filesuffix[],char fileextension[] ) {
	char filename[200];
	strcpy( filename,dir );
	strcat( filename,fileprefix );
	sprintf( filesuffix,"%i",i );
	strcat( filename,filesuffix );
	strcat( filename,fileextension );
	fdetail[i] = fopen( filename,"w+" );
	if( !fdetail[i] ){ // file couldn't be opened
		printf( "Error: File '%s'could not be opened.\n",filename );
		exit( 1 );
	}
	outheader( fdetail[i],i );
	fprintf( fdetail[i],"SPECIES: %i\n",i );
	coordheader( fdetail[i] );
}
void opencoarse( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	coarseheader( *f );
}
void openavvel( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
// 	avvelheader( *favvel );
	avvelWithGradVelheader( *f );
}
void openorder( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	orderheader( *f );
}
void openorderQ( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	orderQheader( *f );
}
void openorderQK( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	orderQKheader( *f );
}
void openavs( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	avsheader( *f );
}
void opendensSTD( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	densheader( *f );
}
void openhistVel( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histVelheader( *f );
}
void openhistSpeed( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histSpeedheader( *f );
}
void openhistVort( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histVortheader( *f );
}
void openhistEnstrophy( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histEnstrheader( *f );
}
void openhistDir( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histDirheader( *f );
}
void openhistS( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histSheader( *f );
}
void openhistDens( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	histNheader( *f );
}
void openavenstrophy( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	avenstrophyheader( *f );
}
void openflow( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	flowheader( *f );
}
void openenergy( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	energyheader( *f );
}
void openenergyfield( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	energyfieldheader( *f );
}
void openenergyneighbours( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	energyneighboursheader( *f );
}
void opensynopsis( FILE **fsynopsis,char dir[],int firsttime ) {
	char filename[200];
	char filesynopsis[]="synopsis";
	char fileextension[]=".dat";

	strcpy( filename,dir );
	strcat( filename,filesynopsis );
	strcat( filename,fileextension );
	if( firsttime ) *fsynopsis = fopen( filename,"w+" );
	else *fsynopsis = fopen( filename,"a" );
	if( !*fsynopsis ){ // file couldn't be opened
		printf("Error: File '%s'could not be opened.\n",filename );
		exit( 1 );
	}
	if( firsttime ) outheader( *fsynopsis,0 );
}
void opentraj( int bc,FILE *fsolids[],char dir[],char filesolids[],char filesuffix[],char fileextension[] ){
	char filename[200];

	strcpy( filename,dir );
	strcat( filename,filesolids );
	sprintf (filesuffix,"%i",bc );
	strcat( filename,filesuffix );
	strcat( filename,fileextension );
	fsolids[bc] = fopen( filename,"w+" );
	if( !fsolids[bc] ){ // file couldn't be opened
		printf( "Error: File '%s'could not be opened.\n",filename );
		exit( 1 );
	}
	outheader( fsolids[bc],bc );
	solidsheader( fsolids[bc] );
}
void openplace( int i,FILE *fin[],char dir[],char fileprefix[],char filesuffix[16],char fileextension[] ){
	char filename[200];

	strcpy( filename,dir );
	strcat( filename,fileprefix );
	sprintf( filesuffix,"%i",i );
	strcat( filename,filesuffix );
	strcat( filename,fileextension );
	fin[i] = fopen( filename,"r" );
	if( !fin[i] ) { // file couldn't be opened
		printf( "Error: File '%s'could not be opened.\n",filename );
		exit( 1 );
	}
}
void opencorr( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	corrheader( *f );
}
void openenergyspect( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	energyspectheader( *f );
}
void openenstrophyspect( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	enstrophyspectheader( *f );
}
void opendefect( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	if(DIM==_3D) printf("Warning: Topological charge field is only outputted for 2D!\n");
	else defectheader( *f );
}
void openmultiphase( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	multiphaseheader( *f );
}
void openpressure( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	pressureheader( *f );
}
void openbinder( FILE **f,char dir[],char fname[],char ext[],int binSize ) {
	openBasic( f,dir,fname,ext );
	binderheader( *f,binSize );
}
void openswimmer( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	swimmerheader( *f );
}
void openswimmerOri( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	swimmeroriheader( *f );
}
void openruntumble( FILE **f,char dir[],char fname[],char ext[] ) {
	openBasic( f,dir,fname,ext );
	runtumbleheader( *f );
}
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ***************** THEORY ***************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

void theory_trans( double *MFP,double *VISC,double *THERMD,double *SDIFF,double *SPEEDOFSOUND,double RA,double FRICCO,double KBT,double dt,double sumM,int RTECH,int SYNOUT,FILE *fsynopsis ) {
/*
	Calculates kinematic viscosity, thermal diffusion coefficient
	and the self-diffusion coefficient. Cell size a =1 is assumed
	This follows:
	Noguchi & Gompper, Transport coefficients of off-lattice mesoscale-hydrodynamics simulation techniques, PRE 78, 016706 (2008)
*/
	double a=1.0;							//MPCD cell size
	double A,B,CM,SM;								//Correlation factors from Table 1 of Nguchi & Gompper
	double VISCKIN,VISCCOL;		//Kinetic and collisional parts of DYNAMIC viscosity
	double THERMDKIN,THERMDCOL;		//Kinetic and collisional parts of thermal diffusion coefficient
	double inv_nDNST;		//Inverse of the number density
	double M;						//Just a value that comes up often
	double inv_DIM;		//Inverse of the dimensionality
	double avMASS = sumM / (double)GPOP;	//Average mass of a particle

	inv_DIM=1.0/(float)DIM;
	inv_nDNST = 1.0/nDNST;
	M = nDNST / (nDNST -1.0 + smrtPow( e,-1.0*nDNST ) );
	//Calculate mean free path
	*MFP = dt * sqrt( KBT/avMASS );
	//Calculate the speed of sound
	*SPEEDOFSOUND = sqrt( 2.0*KBT/avMASS );

	//Calculate the correlation factors
	if(RTECH==ORTHAXIS) {
		A=2.0*inv_DIM*(1.0-cos(RA));
		if(DIM==_2D) B=1.0-cos(2.0*RA);
		else B=(2.0/5.0)*(2.0-cos(RA)-cos(2.0*RA));
	}
	else if(RTECH==ARBAXIS || RTECH==NOHI_ARBAXIS){
		A=2.0*inv_DIM*(1.0-sin(RA)/RA);
		if(DIM==_2D) B=1.0-0.5*sin(2.0*RA)/RA;
		else B=(2.0/5.0)*(2.0-sin(RA)/RA-0.5*sin(RA)/RA);
	}
	else if(RTECH==MPCAT || RTECH==RAT || RTECH==NOHI_MPCAT) {
		A=1.0;
		B=1.0;
	}
	else if(RTECH==LANG || RTECH==RLANG) {
		A=(FRICCO*dt/avMASS)/(1.0+0.5*FRICCO*dt/avMASS);
		B=(2.0*FRICCO*dt/avMASS)/(1.0+0.5*FRICCO*dt/avMASS)/(1.0+0.5*FRICCO*dt/avMASS);
	}
	else {
		// All other versions do not have theoretically known transport coefficients
		A=1.0;
		B=1.0;
	}
	//Calculate the viscosity
	VISCKIN = 0.0;
	VISCCOL = 0.0;
	//Kinetic part of viscosity
	if(RTECH==ORTHAXIS || RTECH==ARBAXIS || RTECH==NOHI_ARBAXIS || RTECH==MPCAT || RTECH==NOHI_MPCAT || RTECH==LANG) {
		//All of the versions without angular-momentum conservation have the same form
		CM=B/M;
		VISCKIN=nDNST*KBT*dt*smrtPow(a,-DIM)*( 1.0/CM-0.5 );
	}
	else if(RTECH==RAT || RTECH==RLANG) {
		//All of the versions with angular-momentum conservation have the same form
		CM=B*(1.0-exp(-nDNST)*(1.0+nDNST)) + (A+inv_DIM*B)*nDNST*exp(-nDNST)/(float)(DIM+2.0);
		CM+=0.5*(A*DIM-0.5*B*(3.0*DIM+2.0))*(1.0-exp(-nDNST)*(1.0+nDNST+0.5*nDNST*nDNST))/nDNST;
		VISCKIN=nDNST*KBT*dt*smrtPow(a,-DIM)*( 1.0/CM-0.5 );
	}
	else {
		//Approximate with only the scaling result
		VISCKIN=nDNST*KBT*dt*smrtPow(a,-DIM);
	}

	//Collisional part of viscosity
	if(RTECH==ORTHAXIS || RTECH==ARBAXIS || RTECH==NOHI_ARBAXIS || RTECH==MPCAT || RTECH==NOHI_MPCAT || RTECH==LANG) {
		//All of the versions without angular-momentum conservation have the same form
		VISCCOL=(A*nDNST/M)*avMASS/(12.0*dt*smrtPow(a,DIM-2));
	}
	else if(RTECH==RAT || RTECH==RLANG) {
		//All of the versions with angular-momentum conservation have the same form
		VISCCOL=(A*avMASS)/(24.0*dt*smrtPow(a,DIM-2))*( nDNST-7.0/5.0+exp(-nDNST)*(7.0/5.0+2.0*nDNST/5.0+(inv_DIM-0.3)*nDNST*nDNST) );
	}
	else {
		//Approximate with only the scaling result
		VISCCOL=avMASS/(dt*smrtPow(a,DIM-2));
	}
	//Total viscosity
	*VISC=VISCKIN+VISCCOL;

	//Calculate the self diffusion coefficient
	if(RTECH==ORTHAXIS || RTECH==ARBAXIS || RTECH==NOHI_ARBAXIS || RTECH==MPCAT || RTECH==NOHI_MPCAT || RTECH==LANG) {
		//All of the versions without angular-momentum conservation have the same form
		SM=A/M;
		*SDIFF = (KBT*dt/avMASS)*(SM-0.5);
	}
	else if(RTECH==RAT || RTECH==RLANG) {
		//All of the versions with angular-momentum conservation have the same form
		SM=0.5*exp(-nDNST)*( 0.5*inv_DIM*nDNST*(DIM-1.0)*(DIM-2.0) + DIM - 1.0 + (DIM+1.0)*inv_nDNST );
		SM+=1.0-0.5*(DIM+1.0)*inv_nDNST;
		SM*=A;
		*SDIFF = (KBT*dt/avMASS)*(SM-0.5);
	}
	else {
		//Approximate with only the scaling result
		*SDIFF = (KBT*dt/avMASS);
	}

	if(RTECH==ORTHAXIS || RTECH==ARBAXIS || RTECH==NOHI_ARBAXIS) {
		//Calculate the thermal diffusion coefficient
		if( DIM ==  _2D ) {
			THERMDKIN = 2. / (1.-cos(RA));
			THERMDKIN -= 1.;
			THERMDKIN += ( 4.*inv_nDNST )*( 1.-0.25 / ( sin( 0.5*RA )*sin( 0.5*RA ) ) );
			THERMDKIN *= 0.5 * KBT * dt;
		}
		else if( DIM == _3D ) {
			THERMDKIN = 0.5 * (2.+cos(RA)) / (1.-cos(RA));
			THERMDKIN += (3.*inv_nDNST) * (0.8-0.25 / (sin( 0.5*RA )*sin( 0.5*RA )));
			THERMDKIN *= KBT*dt;
		}
		else {
			printf( "DIM must be 2 or 3 to calculated thermal diffusivity.\n" );
			THERMDKIN = 0.;
		}
		THERMDCOL = (1./(12.*(DIM+2.)*dt)) * ((nDNST)-0.555555*nDNST*nDNST) * (1.-cos(RA));
		if( THERMDCOL <= 0. ){
			THERMDCOL = 1. + smrtPow( e,-nDNST ) * ( log( nDNST )-1. ) ;
			THERMDCOL -= inv_nDNST;
			THERMDCOL -= smrtPow( inv_nDNST,2. );
			THERMDCOL -= 2. * smrtPow( inv_nDNST,3. );
			THERMDCOL *= (inv_nDNST) / (3.*(DIM+2.)*dt);
			THERMDCOL *= (1.-cos(RA));
		}
	}
	else {
		THERMDKIN = 0.5;
		THERMDCOL = 0.5;
	}
	*THERMD = THERMDKIN+THERMDCOL;
	if( SYNOUT == OUT ) {
		fprintf( fsynopsis,"Mean Free Path: %lf\n",*MFP );
		fprintf( fsynopsis,"Kinematic Viscosity: %lf\n",*VISC );
		fprintf( fsynopsis,"\tkinetic contribution: %lf\n\tcollisional contribution: %lf\n",VISCKIN,VISCCOL );
		fprintf( fsynopsis,"Self Diffusion Coefficient: %lf\n",*SDIFF );
		fprintf( fsynopsis,"Schmidt number: %lf\n",(*VISC)/(*SDIFF) );
		fprintf( fsynopsis,"Speed of sound: %lf\n",*SPEEDOFSOUND );
		fprintf( fsynopsis,"Thermal Diffusion Coefficient: %lf\n",*THERMD );
	}
}

double ndensity( bc WALL[] ) {
/*
   Calculates the number density of the fluid.
   Either per unit volume or area
*/
	int i;
	double V;	//The control volume.
	double D;	//The number density
	V = (double)( XYZ[0] * XYZ[1] * XYZ[2] );
	for( i=0; i<NBC; i++ ) V -= WALL[i].VOL;
	D = GPOP / V;
	return D;
}
double mdensity( bc WALL[],double MASS ) {
/*
   Calculates the mass density of the fluid.
   Either per unit volume or area
*/
	int i;
	double V;	//The control volume.
	double D;	//The mass density
	V = (double)( XYZ[0] * XYZ[1] * XYZ[2] );
	for( i=0; i<NBC; i++ ) V -= WALL[i].VOL;
	D = MASS / V;
	return D;
}
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ***************** ZEROING **************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
void zerovec( double VEC[],int dimension ) {
/*
    Zeros any vector
*/
	int i;
	for( i=0; i<dimension; i++ ) VEC[i]=0.0;
}
void zeroparticles( particleMPC *pp ) {
/*
    Zero everything in all the particles
*/
	int i,d;

	for( i=0; i<GPOP; i++ ) {
		for( d=0; d<_3D; d++ ) {
			(pp+i)->Q[d]=0.0;
			(pp+i)->V[d]=0.0;
			(pp+i)->U[d]=0.0;
			(pp+i)->T[d]=0.0;
		}
	}
}
void zero_bc_var( double *tfrac,double *tdiff,int *g ) {
/*
   This function zeros some of the variables used by BCs
*/
	*tfrac = 0.;
	*tdiff = 0.;
	*g = 0;
}
void zerocnt( double *KBTNOW,double AVNOW[],double *AVS ) {
	int i;
	//Zero counters
	*KBTNOW = 0.;
	*AVS = 0.;
	for( i=0; i<_3D; i++ ) AVNOW[i] = 0.;
}
void zeroHISTVEC( int HIST[_3D][BINS] ) {
	int i,j;
	//Zero counters
	for( j=0; j<BINS; j++ ) for( i=0; i<_3D; i++ ) HIST[i][j] = 0;
}
void zeroHISTSCALAR( int HIST[BINS] ) {
	int j;
	//Zero counters
	for( j=0; j<BINS; j++ ) HIST[j] = 0;
}
void zerocell( cell ***CL ) {
/*
    Zero everything in the cell lists
*/
	int i,j,k,l,m;
	for( i=0; i<XYZ_P1[0]; i++ ) for( j=0; j<XYZ_P1[1]; j++ ) for( k=0; k<XYZ_P1[2]; k++ ) {
		CL[i][j][k].POP = 0;
		CL[i][j][k].MASS = 0.0;
		CL[i][j][k].S = 0.0;
		for( l=0; l<_3D; l++ ) {
			CL[i][j][k].CM[l] = 0.0;
			CL[i][j][k].VCM[l] = 0.0;
			CL[i][j][k].FLOW[l] = 0.0;
			CL[i][j][k].DIR[l] = 0.0;
			for( m=0; m<_3D; m++ ) {
				CL[i][j][k].E[l][m] = 0.0;
				CL[i][j][k].I[l][m] = 0.0;
				CL[i][j][k].Ps[l][m] = 0.0;
				CL[i][j][k].Pc[l][m] = 0.0;
			}
		}
		//The list doesn't have anyone in it yet so it doesn't point anywhere
		CL[i][j][k].pp = NULL;
		CL[i][j][k].MDpp = NULL;
		CL[i][j][k].sp = NULL;
	}
}
void zeroPressureColl( cell *CL ) {
/*
    Zero collisional pressure term
*/
	int l,m;
	for( l=0; l<DIM; l++ ) for( m=0; m<DIM; m++ ) CL->Pc[l][m] = 0.0;
}
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************** INITIALIZING ************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
void initvar( unsigned long *seed,time_t *to,clock_t *co,int *runtime,int *warmtime,double *sumM,double AV[_3D],double avDIR[_3D],spec SP[],double *C,double *S,double RA,double *AVVEL,double KBT,bc WALL[],double MAG[_3D],cell ***CL,particleMPC *pp ) {
/*
   Initializes many program variables and some
   physical parameters.
*/
	int i;

	*seed = RandomSeedSRD (*seed);			//note to tyler

	//init_genrand( time( NULL ) );	//Initialize random number generator.
	*to = time( NULL );
	*co = clock( );
	*runtime = 0;
	*warmtime = 0;
	for( i=0; i<_3D; i++ ) AV[i] = 0.;
	for( i=0; i<_3D; i++ ) avDIR[i] = 0.;
	*sumM = 0.;
	for( i=0; i<NSPECI; i++ ) *sumM += (SP[i].POP) * (SP[i].MASS);
	*C = cos( RA );
	*S = sin( RA );
	*AVVEL = sqrt( DIM*KBT * GPOP / (*sumM) );	//If isotropic velocity dist.
	//Calculate physical parameters of BCs
	for( i=0; i<NBC; i++ ) {		//Set the velocity of the walls to zero
		dim_vol( &WALL[i],XYZ,DIM );
		mominert( &WALL[i],XYZ,DIM );
	}

	nDNST = ndensity( WALL );
	mDNST = mdensity( WALL,*sumM );
	for( i=0; i<DIM; i++ ) MAG[i]/=(nDNST);//Want torque per unit volume so divide field by number density
	//Zero everything in the cell lists
	zerocell( CL );
	//Zero everything in the particles
	zeroparticles( pp );
}

void place( double Q[],int PL,FILE *fin ) {
/*
   This subroutine does the initial placing of the
   particleMPCs. Currently is places them randomly but
   in the future it will have multiple options
*/
	int d;

	if( PL == PRF ) for( d=0; d<DIM; d++ ) Q[d] = genrand_real( ) * XYZ[d];
	else if( PL == READ ) for( d=0; d<DIM; d++ ) {
		if(fscanf( fin, "%lf",&Q[d] ));
		else printf("Warning: Failed to read place input file.\n");
	}
	else{
		printf( "Error: Particle placement distribution unacceptable.\n" );
		exit( 1 );
	}
}
void replace( particleMPC *p ) {
/*
     Randomly place the particleMPC within one SRD cell of its current position
*/
	int i;
	for( i=0; i<DIM; i++ ) p->Q[i] = (double)XYZ[i] * genrand_real( );
}
void push( double V[],double KBT,int PL,double MASS,FILE *fin ) {
/*
   This subroutine does the initial setting of the
   particleMPCs' velocity. Currently is places them
   randomly but in the future it will have multiple options
*/
	int d;
	double normalize;
	if( PL == RANDVEL ) for( d=0; d<DIM; d++ ) V[d] = sqrt( KBT/MASS ) * (2. * genrand_real() - 1.);
							//The particleMPCs have an average speed of sqrt(DIM*KBT/M) BUT each compenent makes this up so v^2=vx^2+vy^2+vz^2 where vx=vy=vz so must divide by sqrt(3) to get component thus no 3.
							//Also for uniform must consider from zero to 2*sqrt(KBT/M).
							//To get negative we make the interval twice as big and subtract off the original max.
	else if(PL == READ) for( d=0; d<DIM; d++ ) {
		if(fscanf( fin, "%lf",&V[d] ));
		else printf("Warning: Failed to read push input file.\n");
	}
	else if( PL == HOMOAXIS ) for( d=0; d<DIM; d++ ) {
		V[d] = sqrt( DIM*KBT/MASS );			//Each particleMPC has the average energy
		V[d] *= genrand_pmOne();			//But could travel in either direction
	}
	else if( PL == FLOCK ) {
		V[0] = sqrt( 2*KBT/MASS );			//Each particleMPC has the average energy
	}
	else if( PL == GAUSS ) for( d=0; d<DIM; d++ ) {
		V[d] = genrand_gaussMB( KBT,MASS );	//Velocity distribution is Gaussian
	}
	else if( PL == HOMO ) {
		for( d=0; d<DIM; d++ ) 	V[d] = genrand_gaussMB( KBT,MASS );	//Velocity distribution is Gaussian (spherically symmetric and separable)
		norm( V,DIM );
		normalize = sqrt((double) DIM *KBT/MASS);
		for( d=0; d<DIM; d++ ) 	V[d] *= normalize;
	}
	else{
		printf( "Error: Particle velocity distribution unacceptable.\n" );
		exit( 1 );
	}
}
void orient( double U[],double Q[],int PL ) {
/*
   This subroutine sets the orientation. Currently is places them
   randomly or all aligned
*/
	int d;
	double noise;
	if( PL==RANDORIENT ) {
		genrand_coneNP( U,pi,DIM );
// 		if(DIM==_2D) U[0]*=genrand_pmOne();
	}
	else if( PL==ALIGNX ) {
		//There is a problem sometimes if it is PERFECTLY aligned so give a bit of noise (+/- 1% max)
		noise = 1.0 + 0.02*( genrand_real()-0.5 );
		U[0]=0.9999*noise;
		U[1]=0.01414178*noise;
		if( DIM>=_3D ) U[2]=0.000008*noise;
	}
	else if( PL==ALIGNY ) {
		//For some reason the collision operation isn't happy if every nematogen is exactly aligned along y
		//so give a bit of noise (+/- 1% max)
		noise = 1.0 + 0.02*( genrand_real()-0.5 );
		U[1]=0.9999*noise;
		U[0]=0.01414178*noise;
		if( DIM>=_3D ) U[2]=0.000008*noise;
	}
	else if( PL==ALIGNZ ) {
		//There is a problem sometimes if it is PERFECTLY aligned so give a bit of noise (+/- 1% max)
		noise = 1.0 + 0.02*( genrand_real()-0.5 );
		U[2]=0.9999*noise;
		U[0]=0.007071*noise;
		U[1]=0.007071*noise;
	}
	else if( PL==ALIGN45 ) {
		if( DIM==_2D ) for( d=0; d<DIM; d++ ) U[d]=1./sqrt(2.);
		else if( DIM==_3D ) for( d=0; d<DIM; d++ ) U[d]=1./sqrt(3.);
		else if( DIM==_1D ) {
			printf( "Error: 1D liquid crystal.\n" );
			exit( 1 );
		}
	}
	else if( PL==PLANEZ ) {
		//Just 2D as far as the orientation is concerned
		genrand_coneNP( U,pi,2 );
	}
	else if( PL==PLANEY ) {
		//Just 2D but switch the ordering of the components to put in XZ-plane
		genrand_coneNP( U,pi,2 );
		U[2]=U[1];
		U[1]=0.0;
	}
	else if( PL==PLANEX ) {
		//Just 2D but switch the ordering of the components to put in YZ-plane
		genrand_coneNP( U,pi,2 );
		U[2]=U[0];
		U[0]=0.0;
	}
	else if( PL==ALIGNTR ) {
		// Orientation points from the origin to the top right corner (furthest point)
		// of the system size, with unit size

		double L2 = 0.0;
		for( d=0; d<DIM; d++ ) L2 += XYZ[d]*XYZ[d];
		L2 = sqrt(L2);

		for( d=0; d<DIM; d++ ) U[d] = XYZ[d]/L2;
	}
	else if( PL==ALIGNDEFECTPAIR ) {
		// Orientation points as if around two oppositely charged defects
		// of the system size, with unit size
		double QP[_2D],QM[_2D];		//Position of the Plus and Minus defects
		double UP[_2D],UM[_2D];		//Orientations of the defects
		double kP=0.5,kM=-0.5;		//Charges of the defects
		double phiP,phiM;			//Director angles at the point Q
		double wP,wM,wY;			//Weightings of U for each defect at the point Q and the walls(linear)
		double r,R;

		//Linear weighting slope
		R=0.5*sqrt( 2.25*((double)XYZ[0])*((double)XYZ[0]) + ((double)XYZ[1])*((double)XYZ[1]) );
		// Do the +1/2 defect on the left side
		QP[0] = ((double)XYZ[0])*0.25;
		QP[1] = ((double)XYZ[1])*0.5;
		phiP = atan2(QP[1]-Q[1],QP[0]-Q[0]);
		UP[0] = cos(kP*phiP);
		UP[1] = sin(kP*phiP);
		r=sqrt( pow(QP[0]-Q[0],2.0)+pow(QP[1]-Q[1],2.0) );
		wP=1.0-r/R;
		// Do the -1/2 defect on the right side
		QM[0] = ((double)XYZ[0])*0.75;
		QM[1] = ((double)XYZ[1])*0.5;
		phiM = atan2(-QM[1]+Q[1],-QM[0]+Q[0]);
		UM[0] = cos(kM*phiM);
		UM[1] = sin(kM*phiM);
		r=sqrt( pow(QM[0]-Q[0],2.0)+pow(QM[1]-Q[1],2.0) );
		wM=1.0-r/R;
		if(Q[1]<=0.5*((double)XYZ[1])) wY=1.0-2.0*Q[1]/(double)XYZ[1];
		else wY=1.0-2.0*((double)XYZ[1]-Q[1])/(double)XYZ[1];

		// Set the director field with a linear weighting from each
		U[0] = wP*UP[0] + wM*UM[0] + wY;
		U[1] = wP*UP[1] + wM*UM[1];
		if( DIM>_2D ) U[2] = 0.0;
		norm( U,DIM );
	}

	else{
		printf( "Error: Particle orientation distribution unacceptable.\n" );
		exit( 1 );
	}
}

int checkplaceMPC( int i,particleMPC *pp,spec SP[],bc WALL[] ) {
	//We must make sure that we check the obstacles that maybe effected by the periodic BC
	double shift[_3D];
	int j,k;

	for( j=0; j<NBC; j++ ) {
		//Zero the shift vector (just in case)
		for( k=0; k<_3D; k++ ) shift[k] = 0.;
		shiftBC( shift,&WALL[j],(pp+i) );
		rotateBC( &WALL[j],(pp+i),0 );
		WALL[j].W = calcW( WALL[j],*(pp+i) );
		rotatebackBC( &WALL[j],(pp+i),0 );
		shiftbackBC( shift,&WALL[j] );
		//If W<=0 then the particleMPC is inside an obstacle and must be replaced
		if( WALL[j].W <= TOL ) {
			replace( (pp+i) );
			//push( (pp+i)->V,KBT,SP[(pp+i)->SPID].VDIST, SP[(pp+i)->SPID].M,NULL );
			i--;
			return i;
		}
	}
	return i;
}
void replacePos_WithCheck( particleMPC *pp,bc WALL[] ) {
	//We must make sure that we check the obstacles that maybe effected by the periodic BC
	double shift[_3D];
	int j,k;
	int flag=1;

	while( flag ) {
		replace( pp );
		//push( pp->V,KBT,SP[pp->SPID].VDIST, SP[pp->SPID].M,NULL );
		flag=0;
		for( j=0; j<NBC; j++ ) {
			//Zero the shift vector (just in case)
			for( k=0; k<_3D; k++ ) shift[k] = 0.;
			shiftBC( shift,&WALL[j],pp );
			rotateBC( &WALL[j],pp,0 );
			WALL[j].W = calcW( WALL[j],*pp );
			rotatebackBC( &WALL[j],pp,0 );
			shiftbackBC( shift,&WALL[j] );
			//If W<=0 then the particleMPC is inside an obstacle and must be replaced
			if( WALL[j].W <= TOL ) flag=1;
		}
	}
}
int checkplace( int IN,particleMPC *pp,spec SP[],bc WALL[],simptr simMD,double KBT,int MDmode ) {
	//We must make sure that we check the obstacles that maybe effected by the periodic BC
	particleMD *atom;
	double d[_3D];
	int i=IN,j;

	i=checkplaceMPC( i,pp,SP,WALL );

	if( MDmode!= noMD ) {
		atom = simMD->atom.items;
		for( j=0; j<simMD->atom.n; j++ ) {
			d[0] = (pp+i)->Q[0] - (atom+j)->rx;
			d[1] = (pp+i)->Q[1] - (atom+j)->ry;
			d[2] = (pp+i)->Q[2] - (atom+j)->rz;

			//Owen made hard code it in instead of call the routine I made. I'm mopey
			if (simMD->pbcond & PBC_COND_x) {
				if (d[0] >= 0.5*XYZ[0]) d[0] -= XYZ[0];
				else if	(d[0] < -0.5*XYZ[0]) d[0] += XYZ[0];
			}
			if (simMD->pbcond & PBC_COND_y) {
				if (d[1] >= 0.5*XYZ[1]) d[1] -= XYZ[1];
				else if	(d[1] < -0.5*XYZ[1]) d[1] += XYZ[1];
			}
			if (simMD->pbcond & PBC_COND_z) {
				if (d[2] >= 0.5*XYZ[2]) d[2] -= XYZ[2];
				else if	(d[2] < -0.5*XYZ[2]) d[2] += XYZ[2];
			}
			if( sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]) <= (simMD->rCut)*1.25 ) {
				replace( (pp+i) );
				//push( (pp+i)->V,KBT,SP[(pp+i)->SPID].VDIST, SP[(pp+i)->SPID].M,NULL );
				i--;
				return i;
			}
		}
	}
	return i;
}

void setcoord( char dir[],spec SP[],particleMPC *pp,double KBT,double AVVEL[],bc WALL[],simptr simMD,int MDmode,int LC ) {
/*
    This subroutine does the initializes the particleMPCs
    coordinates (position, velocity, and orientation. It does so by calling
		place(), push(), and orient()
*/
	int i,j,k=0;
	FILE *fin[NSPECI];
	char fileprefix[] = "placeSP";
	char filesuffix[16] = "0000";
	char fileextension[] = ".inp";

	for( i=0; i<NSPECI; i++ ) for( j=0; j<SP[i].POP; j++ ) {
		(pp+k)->SPID = i;
		k++;
	}
	for( j=0; j<DIM; j++ ) AVVEL[j] = 0.0;
	//Open files incase the particleMPCs' coordinates are being read in
	//Initialize the detailed output files
	for( i=0; i<NSPECI; i++ ) if( SP[i].POP > 0 && SP[i].QDIST == READ ) openplace( i,fin,dir,fileprefix,filesuffix,fileextension );
	//Set  particleMPC coordinates
	for( i=0; i<GPOP; i++ ) {
		//zero (dimensions not used will stay zeroed)
		for( j=0; j<_3D; j++ ) {
			(pp+i)->Q[j] = 0.0;
			(pp+i)->V[j] = 0.0;
			(pp+i)->U[j] = 0.0;
		}
		//Set particleMPC position, velocity and orientation
		place( (pp+i)->Q,SP[(pp+i)->SPID].QDIST, fin[(pp+i)->SPID] );
		//HACK!!!!!!!!!
		//HACK!!!!!!!!!
		//HACK!!!!!!!!!
		//HACK!!!!!!!!!
		//HACK!!!!!!!!!
		//HACK!!!!!!!!!
		//Strip
		// if( SP[(pp+i)->SPID].PHI<0.0 ) (pp+i)->Q[0] *= 0.5;
		// else (pp+i)->Q[0] = 0.5*( (pp+i)->Q[0] + XYZ[0] );
		//Box
		// if( SP[(pp+i)->SPID].PHI<0.0 ) {
		// 	(pp+i)->Q[0] *= 0.25;
		// 	(pp+i)->Q[1] *= 0.25;
		// }
		// else {
		// 	if( (pp+i)->Q[0]<0.25*XYZ[0] && (pp+i)->Q[1]<0.25*XYZ[1] ) {
		// 		(pp+i)->Q[0] += 0.25*XYZ[0];
		// 		(pp+i)->Q[1] += 0.25*XYZ[1];
		// 	}
		// }

		push( (pp+i)->V,KBT,SP[(pp+i)->SPID].VDIST, SP[(pp+i)->SPID].MASS,fin[(pp+i)->SPID] );
		//Shift first mode of the velocity dist by the average velocity (push() centres about zero)
		for( j=0; j<DIM; j++ ) (pp+i)->V[j] += AVVEL[j];

		if( LC>ISOF ) orient( (pp+i)->U,(pp+i)->Q,SP[(pp+i)->SPID].ODIST );
	}
	for( j=0; j<DIM; j++ ) AVVEL[j] /= (double)GPOP;
	//Check particleMPC coordinates
	for( i=0; i<GPOP; i++ ) i = checkplace( i,pp,SP,WALL,simMD,KBT,MDmode );
	//They will all need to stream
	for( i=0; i<GPOP; i++ ) (pp+i)->S_flag = STREAM;

	//Close the input files
	for( i=0; i<NSPECI; i++ ) if( SP[i].POP > 0 && SP[i].QDIST == READ ) fclose( fin[i] );
}
void checkSim( FILE *fsynopsis,int SYNOUT,inputList in,spec *SP,bc *WALL,specSwimmer SS ) {
/*
    This subroutine just checks for odd input
	
*/
	int i,j;

	//Check thermostat
	if( in.TSTECH==BEREND && in.TAU<0.5*in.dt ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Error: TAU = %lf. For Berendsen thermostat TAU must be close to 1*dt.\n",in.TAU );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: TAU = %lf. For Berendsen thermostat TAU must be close to 1*dt.\n",in.TAU );
		exit(1);
	}
	if( in.TSTECH == HEYES && (in.TAU<0.05 || in.TAU>0.3) ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Error: TAU = %lf. For Heyes thermostat must have 0.05<TAU/dt<0.3.\n",in.TAU );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: TAU = %lf. For Heyes thermostat must have 0.05<TAU/dt<0.3.\n",in.TAU );
		exit(1);
	}
	if( ( in.TSTECH==VSC || in.TSTECH==BEREND ) && ( in.RTECH==MPCAT || in.RTECH==RAT ) ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Error: RTECH = %d but TSTECH = %d. For Andersen MPCD, the thermostat must be turned off or apply only to the centre of mass.\n",in.RTECH,in.TSTECH );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: RTECH = %d but TSTECH = %d. For Andersen MPCD, the thermostat must be turned off or apply only to the centre of mass.\n",in.RTECH,in.TSTECH );
		exit(1);
	}
	if( in.TSTECH==HEYES && ( in.RTECH==MPCAT || in.RTECH==RAT ) ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Warning: RTECH = %d but TSTECH = %d. For Andersen MPCD, the including the Heyes thermostat is redundant.\n",in.RTECH,in.TSTECH );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: RTECH = %d but TSTECH = %d. For Andersen MPCD, the including the Heyes thermostat is redundant.\n",in.RTECH,in.TSTECH );
	}
	if( in.RTECH==VICSEK || in.RTECH==CHATE ) for( i=0; i<NSPECI; i++ ) if( SP[i].ACT>1.0 || SP[i].ACT<0.0 ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Error: RTECH = %d (Vicsek or Chate algorithm) but ACT=%lf (the noise) is not within [0,1].\n",in.RTECH,SP[i].ACT );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: RTECH = %d (Vicsek or Chate algorithm) but ACT=%lf (the noise) is not within [0,1].\n",in.RTECH,SP[i].ACT );
		exit(1);
	}
	if( in.TSTECH==MAXV && in.RTECH!=MPCAT ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Warning: RTECH = %d but TSTECH = %d. This thermostat is a maximum velocity controller designed to work with Andersen MPCD.\n",in.RTECH,in.TSTECH );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: RTECH = %d but TSTECH = %d. This thermostat is a maximum velocity controller designed to work with Andersen MPCD.\n",in.RTECH,in.TSTECH );
	}
	// Check number of species
	if( NSPECI > MAXSPECI ) {
		printf("Error: The number of species cannot exceed %d unless definitions.h is altered to increase MAXSPECI.\n", MAXSPECI);
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: The number of species cannot exceed %d unless definitions.h is altered to increase MAXSPECI.\n", MAXSPECI);
		exit(1);
	}
	// Check number of species
	if( NBC > MAXBC ) {
		printf("Error: The number of BCs cannot exceed %d unless definitions.h is altered to increase MAXBC.\n", MAXBC);
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: The number of BCs cannot exceed %d unless definitions.h is altered to increase MAXBC.\n", MAXBC);
		exit(1);
	}
	// Check dimensionality
	if( DIM == _2D ) if ( XYZ[2] != 1 ) {
		printf("Error: In 2D, the simulation must be in the z-plane. Set DZ=1\n");
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: In 2D, the simulation must be in the z-plane. Set DZ=1\n");
		exit(1);
	}
	if( DIM > _3D ) {
		printf("Error: Hyper-dimensional solvent simulations not supported.\n");
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: In 2D, the simulation must be in the z-plane. Set DZ=1\n");
		exit(1);
	}
	if( in.RFRAME == 1 ) for( i=0; i<NSPECI; i++ ) if( SP[i].VDIST == 5 ) {
		printf("Error: Set initialized to FLOCKing but, then subtracted all velocity to zero. Set rFrame to 0\n");
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: Set initialized to FLOCKing but, then subtracted all velocity to zero. Set rFrame to 0\n");
		exit(1);
	}
	// Initialize BC
	//Check that BCs make sense
	for( i=0; i<NBC; i++ ) {
		if( WALL[i].DSPLC==0 ) for( j=0; j<DIM; j++ ) {
			if( fneq(WALL[i].V[j],0.0) ) {
				printf( "Warning:\tBC %d immobile; setting velocity to zero.\n",i );
				if(SYNOUT == OUT) fprintf(fsynopsis,"Warning:\tBC %d immobile; setting velocity to zero.\n",i );
				WALL[i].V[j]=0.0;
			}
			if( fneq(WALL[i].L[j],0.0) ) {
				printf( "Warning:\tBC %d immobile; setting angular velocity to zero.\n",i );
				if(SYNOUT == OUT) fprintf(fsynopsis,"Warning:\tBC %d immobile; setting angular velocity to zero.\n",i );
				WALL[i].L[j]=0.0;
			}
			if( fneq(WALL[i].G[j],0.0) ) {
				printf( "Warning:\tBC %d immobile; setting acceleration to zero.\n",i );
				if(SYNOUT == OUT) fprintf(fsynopsis,"Warning:\tBC %d immobile; setting acceleration to zero.\n",i );
				WALL[i].V[j]=0.0;
			}
		}
		if( WALL[i].ABS==0 ) {
			// The user should likely use the abs() operator for BCs that aren't a plane or don't have even powers
			for( j=0; j<DIM; j++ ) {
				if( fneq(WALL[i].P[j],1.0) ) if( fneq( fmod(WALL[i].P[j],2.0),0.0 ) ) {
					printf( "Warning:\tBC %d does not use abs() but likely should when P[%d]=%lf is not even or unity.\n",i,j,WALL[i].P[j] );
					printf( "\tIf the user means to create a closed obstacle or particle this is likely to generate errors.\n" );
					printf( "\t***However***, this is left to the user's discretion.\n" );
					if(SYNOUT == OUT) {
						fprintf(fsynopsis,"Warning:\tBC %d does not use abs() but likely should when P[%d]=%lf is not even or unity.\n",i,j,WALL[i].P[j] );
						fprintf( fsynopsis,"\tIf the user means to create a closed obstacle or particle this is likely to generate errors.\n" );
						fprintf( fsynopsis,"\t***However***, this is left to the user's discretion.\n" );
					}
				}
			}
		}
	}
	if( !( in.LC==ISOF || in.LC==LCL || in.LC==LCG) ){
		printf( "Error: Unrecognized value of LC=%d.\n",in.LC );
		exit( 1 );
	}
	if( in.LC==LCG ) for( i=0; i<NSPECI; i++ ) if( SP[i].ODIST==RANDORIENT ) {
		printf( "Warning: Using global S (LC=%d) but initiated in isotropic phase. Simulation may not reach nematic phase.\n",in.LC );
		if(SYNOUT == OUT) fprintf( fsynopsis,"Warning: Using global S (LC=%d) but initiated in isotropic phase. Simulation may not reach nematic phase.\n",in.LC );
	}
	//Check that nematogens have non-zero friction
	if( in.LC>ISOF ) for( i=0; i<NSPECI; i++ ) if( feq(SP[i].RFC,0.0) ) {
		printf( "Warning:\tSpecies %d has zero rotational friction coefficient\n",i );
		printf( "\t\tThis is acceptable iff no magnetic field is applied\n");
		if(SYNOUT == OUT) fprintf( fsynopsis,"Warning:\tSpecies %d has zero rotational friction coefficient\n\t\tThis is acceptable iff no magnetic field is applied\n",i );
	}
	//Check that if running as a LC that the LC mean field potential MFPOT is greater than zero
	if( in.LC>ISOF && in.MFPOT<=0.0 ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Error: Running as nematic liquid crystal (LC=%d) but mean field = %lf. Must be greater than 0 or run as isotropic\n",in.LC,in.MFPOT );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: Running as nematic liquid crystal (LC=%d) but mean field = %lf. Must be greater than 0 or run as isotropic\n",in.LC,in.MFPOT );
		exit(1);
	}
// 	if( in.LC==LCG || in.LC==LCL ) for( i=0; i<NSPECI; i++ ) if( SP[i].RFC*SP[i].TUMBLE>1.0 ) {
// 		printf( "Warning: Tumbling parameter (%lf) times rotational friction coefficient (%lf) greater than 1. Simulation likely to crash.\n",SP[i].TUMBLE,SP[i].RFC );
// 		if(SYNOUT == OUT) fprintf( fsynopsis,"Warning: Tumbling parameter (%lf) times rotational diffusion coefficient (%lf) greater than 1. Simulation likely to crash.\n",SP[i].TUMBLE,SP[i].RFC );
// 	}
	if( in.LC!=ISOF ) for( i=0; i<NSPECI; i++ ) if( SP[i].ODIST==ALIGNZ && DIM<_3D ) {
		#ifdef DBG
			if( DBUG >= DBGWARN ) printf( "Error: \tSpecies %d orientations initialized along z-axis but not 3D\n",i );
		#endif
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: \tSpecies %d orientations initialized along z-axis but not 3D\n",i );
		exit(1);
	}
	if( in.LC!=ISOF ) if( in.RTECH==CHATE || in.RTECH==CHATE_MPCAT || in.RTECH==CHATE_LANG || in.RTECH==DIPOLE_VCM ) {
		printf( "Error: RTECH is a Chate-style active nematic but LC should be %d (LC=%d and RTECH=%d).\n",ISOF,in.LC,in.RTECH );
		if(SYNOUT == OUT) fprintf( fsynopsis,"Error: RTECH is a Chate-style active nematic but LC should be %d (LC=%d and RTECH=%d).\n",ISOF,in.LC,in.RTECH );
		exit( 1 );
	}
	if( (in.RTECH==DIPOLE_DIR_SUM || in.RTECH==DIPOLE_DIR_AV) && in.LC==ISOF ) {
		printf( "Error: RTECH is a dipole force along the director (RTECH=%d) but LC is for an isotropic fluid (LC=%d).\n",in.RTECH,in.LC );
		if(SYNOUT == OUT) fprintf( fsynopsis, "Error: RTECH is a dipole force along the director (RTECH=%d) but LC is for an isotropic fluid (LC=%d).\n",in.RTECH,in.LC );
		exit( 1 );
	}
	// Check binary fluid
	// if( in.binaryDELTA > 1.0 || in.binaryDELTA < 0.0 ) {
	// 	printf("Warning: The binary fluid fractional energy difference (binaryDELTA) should be between zero and unity.\n");
	// 	if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: The binary fluid fractional energy difference (binaryDELTA) should be between zero and unity.\n");
	// 	exit(1);
	// }
	// Check swimmers
	if( SS.sizeShrink>1.0 ) {
		printf("Warning: The swimmers grow (smaller rotational diffusion) during tumble phase rather than shrink (sizeShrink=%lf).\n",SS.sizeShrink);
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: The swimmers grow (smaller rotational diffusion) during tumble phase rather than shrink (sizeShrink=%lf).\n",SS.sizeShrink);
	}
	if( feq(SS.DS,0.0) ) {
		printf("Warning: The swimmers' dipole strength is zero so they do not geneate dipolar flow.\n");
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: The swimmers' dipole strength is zero so they do not geneate dipolar flow.\n");
	}
	if( !feq(fabs(SS.TS),0.0) && DIM<_3D ) {
		printf("Warning: The swimmers' rotlet dipole strength is non-zero but the simulation is 2D so no rotlet possible.\n");
		if(SYNOUT == OUT) fprintf(fsynopsis,"Warning: The swimmers' rotlet dipole strength is non-zero but the simulation is 2D so no rotlet possible.\n");
	}
	if( SS.runTime<0.0 || SS.tumbleTime<0.0 ) {
		printf("Error: The run or tumble time must be greater than zero (or zero for no run/tumble dynamics).\n");
		if(SYNOUT == OUT) fprintf(fsynopsis,"Error: The run or tumble time must be greater than zero (or zero for no run/tumble dynamics).\n");
		exit( 1 );
	}
}

void initOutput( char op[],outputFlagsList *outFlag,outputFilesList *outFile,inputList in,spec *SP, bc WALL[] ) {

	int i;
	char filecoarse[]="coarsegrain";
	char fileavvel[]="avVel";
	char fileorder[]="directorfield";
	char fileorderQ[]="ordertensor";
	char fileorderQK[]="recipOrder";
	char fileavs[]="avS";
	char filedensSTD[]="densSTD";
	char fileenstrophy[]="avEnstrophy";
	char fileflow[]="flowfield";
	char filesolids[]="solidtraj";
	char filedefect[]="topochargefield";
	char fileprefix[]="detailedSP";
	char filehistVel[]="distVel";
	char filehistVort[]="distVort";
	char filehistDens[]="distDens";
	char filehistS[]="distS";
	char filehistDir[]="distDir";
	char filehistSpeed[]="distSpeed";
	char filehistEnstr[]="distEnstrophy";
	char fileenergy[]="energy";
	char fileenergyfield[]="enfield";
	char fileenneighbours[]="enneighbours";
	char filecorrVV[]="corrVelVel";
	char filecorrNN[]="corrDirDir";
	char filecorrDD[]="corrDensDens";
	char filecorrSS[]="corrOrderOrder";
	char filecorrPP[]="corrPhiPhi";
	char filecorrWW[]="corrVortVort";
	char fileenergyspect[]="energySpectrum";
	char fileenstrophyspect[]="enstrophySpectrum";
	char fileBinder[]="binderCumulant";
	char fileswimmer[]="swimmers";
	char fileswimmerori[]="swimmersOri";
	char fileruntumble[]="runtumble";
	char filemultiphase[]="multiphase";
	char filepressure[]="pressure";

	char filesuffix[]="0000";
	char fileextension[]=".dat";

	#ifdef DBG
		if( DBUG >= DBGINIT ) printf("Initialization\n");
		if( DBUG >= DBGINIT ) printf("Initialize Output Files\n");
	#endif
	//Don't bother with LC stuff if it's not being used
	if(in.LC==ISOF) {
		outFlag->ORDEROUT=0;
		outFlag->QTENSOUT=0;
		outFlag->QKOUT=0;
		outFlag->CNNOUT=0;
		outFlag->CSSOUT=0;
	}
	//Initialize the detailed output files
	if( (outFlag->TRAJOUT)>=OUT ) for(i=0;i<NSPECI;i++) if(SP[i].POP>=1) opendetails( i,outFile->fdetail,op,fileprefix,filesuffix,fileextension );
	//Initialize the course grained output file
	if( (outFlag->COAROUT)>=OUT ) opencoarse( &(outFile->fcoarse),op,filecoarse,fileextension );
	//Initialize the global average velocity MPCD output file
	if( (outFlag->AVVELOUT)>=OUT ) openavvel( &(outFile->favvel),op,fileavvel,fileextension );
	//Initialize the director output file
	if( (outFlag->ORDEROUT)>=OUT ) openorder( &(outFile->forder),op,fileorder,fileextension );
	//Initialize the tensor order parameter output file
	if( (outFlag->QTENSOUT)>=OUT ) openorderQ( &(outFile->forderQ),op,fileorderQ,fileextension );
	if( (outFlag->QKOUT)>=OUT ) openorderQK( &(outFile->forderQK),op,fileorderQK,fileextension );
	//Initialize the scalar order parameter output file
	if( (outFlag->AVSOUT)>=OUT ) openavs( &(outFile->favs),op,fileavs,fileextension );
	//Initialize the scalar order parameter output file
	if( (outFlag->DENSOUT)>=OUT ) opendensSTD( &(outFile->fdensSTD),op,filedensSTD,fileextension );
	//Initialize the enstrophy output file
	if( (outFlag->ENSTROPHYOUT)>=OUT ) openavenstrophy( &(outFile->fenstrophy),op,fileenstrophy,fileextension );
	//Initialize the flow field output file
	if( (outFlag->FLOWOUT)>=OUT ) openflow( &(outFile->fflow),op,fileflow,fileextension );
	//Initialize the distribution output files
	if( (outFlag->HISTVELOUT)>=OUT ) openhistVel( &(outFile->fhistVel),op,filehistVel,fileextension );
	if( (outFlag->HISTSPEEDOUT)>=OUT ) openhistSpeed( &(outFile->fhistSpeed),op,filehistSpeed,fileextension );
	if( (outFlag->HISTVORTOUT)>=OUT ) openhistVort( &(outFile->fhistVort),op,filehistVort,fileextension );
	if( (outFlag->HISTENSTROUT)>=OUT ) openhistEnstrophy( &(outFile->fhistEnstr),op,filehistEnstr,fileextension );
	if( (outFlag->HISTDIROUT)>=OUT ) openhistDir( &(outFile->fhistDir),op,filehistDir,fileextension );
	if( (outFlag->HISTSOUT)>=OUT ) openhistS( &(outFile->fhistS),op,filehistS,fileextension );
	if( (outFlag->HISTNOUT)>=OUT ) openhistDens( &(outFile->fhistDens),op,filehistDens,fileextension );
	//Initialize the energy output files
	if( (outFlag->ENOUT)>=OUT ) openenergy( &(outFile->fenergy),op,fileenergy,fileextension );
	//Initialize the energy field output files
	if( (outFlag->ENFIELDOUT)>=OUT ) openenergyfield( &(outFile->fenergyfield),op,fileenergyfield,fileextension );
	//Initialize the energy from neighbours output files
	if( (outFlag->ENNEIGHBOURS)>=OUT ) openenergyneighbours( &(outFile->fenneighbours),op,fileenneighbours,fileextension );
	//Initialize the spatial correlation output files
	if( (outFlag->CVVOUT)>=OUT ) opencorr( &(outFile->fcorrVV),op,filecorrVV,fileextension );
	if( (outFlag->CWWOUT)>=OUT ) opencorr( &(outFile->fcorrWW),op,filecorrWW,fileextension );
	if( (outFlag->CDDOUT)>=OUT ) opencorr( &(outFile->fcorrDD),op,filecorrDD,fileextension );
	if( (outFlag->CPPOUT)>=OUT ) opencorr( &(outFile->fcorrPP),op,filecorrPP,fileextension );
	if( (outFlag->CNNOUT)>=OUT ) opencorr( &(outFile->fcorrNN),op,filecorrNN,fileextension );
	if( (outFlag->CSSOUT)>=OUT ) opencorr( &(outFile->fcorrSS),op,filecorrSS,fileextension );
	//Initialize the spectrum output files
	if( (outFlag->ENERGYSPECTOUT)>=OUT ) openenergyspect( &(outFile->fenergyspect),op,fileenergyspect,fileextension );
	if( (outFlag->ENSTROPHYSPECTOUT)>=OUT ) openenstrophyspect( &(outFile->fenstrophyspect),op,fileenstrophyspect,fileextension );
	//Initialize the Binder cumulant output files
	if( (outFlag->BINDER)>=OUT ) openbinder( &(outFile->fbinder),op,fileBinder,fileextension,outFlag->BINDERBIN );
	//Initialize swimmer output files
	if( (outFlag->SWOUT)>=OUT ) openswimmer( &(outFile->fswimmers),op,fileswimmer,fileextension );
	//Initialize swimmer output files
	if( (outFlag->SWORIOUT)>=OUT ) openswimmerOri( &(outFile->fswimmersOri),op,fileswimmerori,fileextension );
	//Initialize swimmer run/tumble output files
	if( (outFlag->RTOUT)>=OUT ) openruntumble( &(outFile->fruntumble),op,fileruntumble,fileextension );
	//Initialize the synopsis output files
	if( (outFlag->SYNOUT)>=OUT ) opensynopsis( &(outFile->fsynopsis),op,1 );
	//Initialize the solids' trajectories (or BC motion) output files
	if( (outFlag->SOLOUT)>=OUT ) for ( i=0; i<NBC; i++ ) if ( WALL[i].DSPLC ) opentraj(i,outFile->fsolids,op,filesolids,filesuffix,fileextension);
	//Initialize the defect trajectories output file
	if( (outFlag->DEFECTOUT)>=OUT ) opendefect( &(outFile->fdefects),op,filedefect,fileextension );
	//Initialize the phi/color/species-type field output file
	if( (outFlag->SPOUT)>=OUT ) openmultiphase( &(outFile->fmultiphase),op,filemultiphase,fileextension );
	//Initialize the pressure field output file
	if( (outFlag->PRESOUT)>=OUT ) openpressure( &(outFile->fpressure),op,filepressure,fileextension );

	if( (outFlag->SYNOUT)==OUT) {
		fprintf(outFile->fsynopsis,"Output:\n" );
		fprintf(outFile->fsynopsis,"\tTrajectories:\t%d\n",outFlag->TRAJOUT);
		fprintf(outFile->fsynopsis,"\tCoarse-Grained Flow:\t%d\n",outFlag->COAROUT);
		fprintf(outFile->fsynopsis,"\tGlobal Average velocity:\t%d\n",outFlag->AVVELOUT);
		fprintf(outFile->fsynopsis,"\tFlow:\t\t%d\n",outFlag->FLOWOUT);
		fprintf(outFile->fsynopsis,"\tPrint distributions:\n");
		fprintf(outFile->fsynopsis,"\t\tVel: %d\n\t\tSpeed: %d\n\t\tVorticity: %d\n\t\tEnstrophy: %d\n\t\tDirector: %d\n\t\tScalar order parameter: %d\n\t\tDensity: %d\n",outFlag->HISTVELOUT,outFlag->HISTSPEEDOUT,outFlag->HISTVORTOUT,outFlag->HISTENSTROUT,outFlag->HISTDIROUT,outFlag->HISTSOUT,outFlag->HISTNOUT);
		fprintf(outFile->fsynopsis,"\tGlobal average scalar order parameter:\t%d\n",outFlag->AVSOUT);
		fprintf(outFile->fsynopsis,"\tAverage enstrophy:\t%d\n",outFlag->ENSTROPHYOUT);
		fprintf(outFile->fsynopsis,"\tStandard deviation of density:\t%d\n",outFlag->DENSOUT);
		fprintf(outFile->fsynopsis,"\tDirector and scalar order parameter:\t%d\n",outFlag->ORDEROUT);
		fprintf(outFile->fsynopsis,"\tOrder parameter tensor:\t%d\n",outFlag->QTENSOUT);
		fprintf(outFile->fsynopsis,"\tOrder parameter tensor (reciprocal space):\t%d\n",outFlag->QKOUT);
		fprintf(outFile->fsynopsis,"\tEnergy:\t\t%d\n",outFlag->ENOUT);
		fprintf(outFile->fsynopsis,"\tEnergy field:\t\t%d\n",outFlag->ENFIELDOUT);
		fprintf(outFile->fsynopsis,"\tEnergy neighbours:\t\t%d\n",outFlag->ENNEIGHBOURS);
		if(DIM==_3D) fprintf(outFile->fsynopsis,"\tTopological charge field:\t\t%d\tNot outputted in 3D!\n",outFlag->DEFECTOUT);
		else fprintf(outFile->fsynopsis,"\tTopological charge field:\t\t%d\n",outFlag->DEFECTOUT);
		fprintf(outFile->fsynopsis,"\tPhi/colour/species-type field:\t\t%d\n",outFlag->SPOUT);
		fprintf(outFile->fsynopsis,"\tPressure field:\t\t%d\n",outFlag->PRESOUT);
		fprintf(outFile->fsynopsis,"\tVelocity-velocity correlation:\t\t%d\n",outFlag->CVVOUT);
		fprintf(outFile->fsynopsis,"\tDirector-director correlation:\t\t%d\n",outFlag->CNNOUT);
		fprintf(outFile->fsynopsis,"\tVorticity-vorticity correlation:\t\t%d\n",outFlag->CWWOUT);
		fprintf(outFile->fsynopsis,"\tDensity-density correlation:\t\t%d\n",outFlag->CDDOUT);
		fprintf(outFile->fsynopsis,"\tOrder-order correlation:\t\t%d\n",outFlag->CSSOUT);
		fprintf(outFile->fsynopsis,"\tPhase-phase (binary fluid) correlation:\t\t%d\n",outFlag->CPPOUT);
		fprintf(outFile->fsynopsis,"\tBinder cumulant:\t\t%d --- bin size:\t\t%d\n",outFlag->BINDER,outFlag->BINDERBIN);
		fprintf(outFile->fsynopsis,"\tSolid BC:\t%d\n",outFlag->SOLOUT);
		fprintf(outFile->fsynopsis,"\tSwimmers:\t%d\n",outFlag->SWOUT);
		fprintf(outFile->fsynopsis,"Files initialized.\n" );
	}
}

void initializeSIM( cell ***CL,particleMPC *SRDparticles,spec SP[],bc WALL[],simptr simMD,specSwimmer *specS,swimmer *swimmers,int argc, char* argv[],inputList *in,time_t *to,clock_t *co,int *runtime,int *warmtime,double *AVVEL,kinTheory *theory,double *KBTNOW,double *AVS,double *S4,double *stdN,double AVNOW[_3D],double AVV[_3D],double avDIR[_3D], outputFlagsList outFlags,int MDmode,FILE *fsynopsis,char ip[] ) {
/*
   Initializes simulation
*/
	int i,j;
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf("\tInitialize Parameters\n");
	#endif
	initvar( &(in->seed),to,co,runtime,warmtime,&(theory->sumM),AVNOW,avDIR,SP,&(in->C),&(in->S),in->RA,AVVEL,in->KBT,WALL,in->MAG,CL,SRDparticles );

	maxXYZ=(int) sqrt( (double)(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]+XYZ[2]*XYZ[2]) );

	in->GRAV_FLAG = 0;
	if( fneq(in->GRAV[0],0.0) || fneq(in->GRAV[1],0.0) || fneq(in->GRAV[2],0.0) ) in->GRAV_FLAG = 1;
	for( i=0; i<NBC; i++ ) if( fneq((WALL+i)->G[0],0.0) || fneq((WALL+i)->G[1],0.0) || fneq((WALL+i)->G[2],0.0) ) in->GRAV_FLAG = 1;
	// Flag whether or not to do BC re-orientations
	for( i=0; i<NBC; i++ ) {
		(WALL+i)->REORIENT = 1;
		// If the colloid cannot move and doesn't have an initial non-zero orientation, then never bother with rotations
		if( (WALL+i)->DSPLC==0 && feq((WALL+i)->O[0],0.0) && feq((WALL+i)->O[1],0.0) && feq((WALL+i)->O[2],0.0) ) (WALL+i)->REORIENT = 0;
		// If the colloid has circular symmetry, then never bother with rotations
		else if( DIM==_2D && feq((WALL+i)->P[0],2.0) && feq((WALL+i)->P[1],2.0) && feq((WALL+i)->A[0],(WALL+i)->A[1]) ) (WALL+i)->REORIENT = 0;
		// If the colloid has spherical symmetry, then never bother with rotations
		else if( feq((WALL+i)->P[0],2.0) && feq((WALL+i)->P[1],2.0) && feq((WALL+i)->P[2],2.0) && feq((WALL+i)->A[0],(WALL+i)->A[1]) && feq((WALL+i)->A[1],(WALL+i)->A[2]) ) (WALL+i)->REORIENT = 0;
	}
	in->MAG_FLAG = 0;
	if( fneq(in->MAG[0],0.0) || fneq(in->MAG[1],0.0) || fneq(in->MAG[2],0.0) ) in->MAG_FLAG = 1;

	if( in->TSTECH==MAXV ) for( j=0; j<DIM; j++ ) {
		AVV[j]=in->GRAV[j];
		in->GRAV[j]=0.0;
	}
	else for( j=0; j<DIM; j++ ) AVV[j]=0.0;

	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "\tInitialize System\n" );
		if( DBUG >= DBGINIT ) printf( "\tInitialize MPC\n" );
	#endif
	//Intialize positions, velocity and angular velocity
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "\tPlace MPC particle\n" );
	#endif
	setcoord( ip,SP,SRDparticles,in->KBT,AVV,WALL,simMD,MDmode,in->LC );
	//Intialize positions and orientations of swimmers
	if( NS>0 ) {
		#ifdef DBG
			if( DBUG >= DBGINIT ) printf( "\tPlace swimmers\n" );
		#endif
		setswimmers( specS,swimmers,WALL,in->stepsMD,in->dt );
	}

	if(outFlags.SYNOUT == OUT) fprintf(fsynopsis,"\nMPC particles placed.\n" );
	//Calculate the theoretical properties of the SRD gas
	theory_trans( &(theory->MFP),&(theory->VISC),&(theory->THERMD),&(theory->SDIFF),&(theory->SPEEDOFSOUND),in->RA,in->FRICCO,in->KBT,in->dt,theory->sumM,in->RTECH,outFlags.SYNOUT,fsynopsis );
		//Zero counters
	zerocnt( KBTNOW,AVNOW,AVS );
	if( in->TSTECH==MAXV ) for( j=0; j<DIM; j++ ) AVNOW[j]=AVV[j];
	//Bin particles
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "\tBin MPC particles\n" );
	#endif
	binin( SRDparticles,CL );
	bininSwimmers( *specS,swimmers,CL );
	if(outFlags.SYNOUT == OUT) fprintf(fsynopsis,"\nMPC particles binned for first time.\n" );
	if( MDmode ) bininMD( simMD,CL );
	localPROP( CL,SP,*specS,in->RTECH,in->LC );
	*S4=0.;
	*stdN=0.;
	if( outFlags.AVSOUT>=OUT ) {
		*AVS = avOrderParam( SRDparticles,in->LC,avDIR );
		*S4 = avS4( SRDparticles,in->LC,avDIR );
		*stdN = stdNum( CL,GPOP,XYZ,XYZ_P1 );
		}
	/* ****************************************** */
	/* ******** GALILEAN TRANSFORMATION ********* */
	/* ****************************************** */
	*KBTNOW = TEMP( SRDparticles,SP,WALL,AVV );
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "Actual System Temperature before moving to rest frame: %lf\n",*KBTNOW );
	#endif
	if( in->RFRAME==1 && GPOP<=2 ) {
		in->RFRAME=0;
		printf( "Simulation of 2 particles or less in the rest frame is boring\nDon't do Galilean Transformation\n" );
	}
	
	#ifdef DBG
		if( DBUG >= DBGINIT ) {
			printf( "Inputted System Temperature (units of KB): %lf\n",in->KBT );
			if( in->RFRAME==1) printf( "Galilean Transformation to System Rest Frame\n" );
		}
	#endif
	//Do the Galilean transformation of the system to it's rest frame i.e. remove system's net momentum
	if( in->RFRAME ) galileantrans( SRDparticles,WALL,simMD,SP,in->KBT,AVV,GPOP,NBC,MDmode,DIM );
	//Now that the initial shift is done we use RFRAME to signal when it should happen periodically (with zeroNetMom)
	//But don't want to do it if accelerating, duh
	if( !(in->zeroNetMom) ) in->RFRAME = 0;
	else if( in->RFRAME ) if( fneq(in->GRAV[0],0.0) || fneq(in->GRAV[1],0.0) || fneq(in->GRAV[2],0.0) ) in->RFRAME = 0;
	if(outFlags.SYNOUT == OUT) fprintf(fsynopsis,"\nGalilean transformation completed.\n" );
	/* ****************************************** */
	/* *********** TEMPERATURE SCALING ********** */
	/* ****************************************** */
	*KBTNOW = TEMP( SRDparticles,SP,WALL,AVV );
	#ifdef DBG
		if( DBUG >= DBGINIT ) {
			printf( "System Temperature before scaling: %lf\n",*KBTNOW );
			printf( "Apply thermostat %d.\n",in->TSTECH );
		}
	#endif
	//Scale temperature TSTECH=VSC here so that temperature moved instantaneously to initialized value
	scaleT( in->KBT,*KBTNOW,in->dt,in->TAU,AVV,AVNOW,VSC,SP,in->LC,WALL,SRDparticles,CL );
	*KBTNOW = TEMP( SRDparticles,SP,WALL,AVV );
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "System Temperature after moving to rest frame and scaling: %lf\n",*KBTNOW );
	#endif
	if(outFlags.SYNOUT == OUT) fprintf(fsynopsis,"Temperature rescaled to input value.\n" );
	avVel( CL,AVNOW );
}

void initializeRecovery( cell ***CL,particleMPC *SRDparticles,spec SP[],specSwimmer specS, int RTECH,int LC,int MDmode,int SYNOUT,FILE *fsynopsis ) {
	//int i;
	if(SYNOUT == OUT) fprintf(fsynopsis,"\nSimulation recovered from checkpoint.\n" );
	// MD isn't checkpointed so can't recover MD simulation
	if( MDmode != noMD ) {
		printf("Error: Cannot recover MD simulation.\n");
		exit(EXIT_FAILURE);
	}
	//maxXYZ=0;
	//for( i=0;i<_3D;i++ ) if( XYZ[i]>maxXYZ ) maxXYZ=XYZ[i];
	maxXYZ=(int) sqrt( (double)(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]+XYZ[2]*XYZ[2]) );
	//Bin particles
	#ifdef DBG
		if( DBUG >= DBGINIT ) printf( "\tBin MPC particles\n" );
	#endif
	binin( SRDparticles,CL );
	if(SYNOUT == OUT) fprintf(fsynopsis,"\nMPC particles binned for first time.\n" );
	localPROP( CL,SP,specS,RTECH,LC );
}
