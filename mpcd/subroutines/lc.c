///
/// @file
/// @brief Functions needed to include liquid crystal in the MPCD code (either as solvent particles, or individual MD rod-like particles).
///
///

# include<math.h>
# include<stdio.h>
# include<stdlib.h>

# include "../headers/definitions.h"
# include "../headers/globals.h"
# include "../headers/SRDclss.h"
# include "../headers/rand.h"
# include "../headers/pout.h"
# include "../headers/mtools.h"
# include "../headers/bc.h"
# include "../headers/lc.h"
# include "../headers/mpc.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* *************** COLLISIONS *************** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
///
/// @brief Generate a random, normalized vector by Maier-Saupe distribution.
///
/// In this routine, different cases for different dimensions (1D, 2D and 3D) and posible limit cases are identified.
/// Based on that, respective subroutines are used to generate a vector.
/// Depending on whether BUS=effM*S/KBT is small or large, different functions are called in the 3-D case:
/// - When BUS is large then a gaussian approximation for the angle works in 3D (2D is turned off currently);
/// - When BUS is small an expansion of the exponent works in 3D (BUT NOT 2D, unfortunately).
/// All other cases use the Metropolis algorithm, which is also always used in the 2-D and 1-D case.
/// For 1D, the orientation MUST lie on the axis regardless of anything else.
///
/// @param DIR The director of the local cell.
/// @param rotAx The axis that the x-axis must be rotated about compared the local director.
/// @param rotAngle The angle between the x-axis and the local director.
/// @param U The vector components generated with this function.
/// @param KBT The temperature.
/// @param S  The nematic magnitude of the local cell.
/// @param effM The molecular field potential (scaled with dimension).
/// @see genrand_maierSaupeGAUSS_3D()
/// @see genrand_maierSaupeEXP_3D()
/// @see genrand_maierSaupeMetropolis_3D()
/// @see genrand_maierSaupeMetropolis_2D()
///
void genrand_maierSaupe( double DIR[],double rotAx[],double rotAngle,double U[],double KBT,double S,double effM ) {
	double BUS=effM*S/KBT;
	if(DIM==_3D) {
		if( BUS>BUSMAX ) genrand_maierSaupeGAUSS_3D( rotAx,rotAngle,U,KBT,S,effM );
		else if( BUS<BUSMIN ) genrand_maierSaupeEXP_3D( rotAx,rotAngle,U,KBT,S,effM );
		else genrand_maierSaupeMetropolis_3D( DIR,rotAx,rotAngle,U,KBT,S,effM );
	}
	else if(DIM==_2D) {
		genrand_maierSaupeMetropolis_2D( DIR,rotAx,rotAngle,U,KBT,S,effM );
		// if( BUS>BUSMAX ) genrand_maierSaupeGAUSS_2D( rotAx,rotAngle,U,KBT,S,effM );
// 		//else if( BUS<BUSMIN ) genrand_maierSaupeEXP_2D( DIR,rotAx,rotAngle,U,KBT,S,effM );
// 		else genrand_maierSaupeMetropolis_2D( DIR,rotAx,rotAngle,U,KBT,S,effM );
	}
	else if(DIM==_1D) {
		// In 1D, the orientation MUST lie on the axis regardless of anything else
		U[0]=genrand_pmOne();
	}
	else printf("Warning: genrand_maierSaupe() only programmed for DIM={3,2}, not DIM=%d\n",DIM);
}

///
/// @brief Generate a random, normalized vector by Maier-Saupe distribution in the large effM*S/KBT limit and 3D.
///
/// The function generates a random axis, and performs the Rodrigues rotation algorithm on rotAx using rotAngle. 
/// The result is stored in U[].
///
/// @param rotAx The axis that will be rotated.
/// @param rotAngle The angle of the rotation.
/// @param U The vector components generated with this function.
/// @param KBT The temperature.
/// @param S  The nematic magnitude of the local cell.
/// @param effM The molecular field potential (scaled with dimension).
///
void genrand_maierSaupeGAUSS_3D( double rotAx[],double rotAngle,double U[],double KBT,double S,double effM ) {
	double std,phi,theta,uz;
	double newU[_3D]={0.0};
	int d;

	for( d=0; d<DIM; d++ ) newU[d]=U[d];
	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Old orientation (gauss): " );
			pvec( U,_3D );
		}
	#endif

	//Generate a new angle
	std=1./sqrt(2.*effM*S/KBT);
	theta=genrand_gaussGen(0.0,std);
// 	uz=cos(theta);
	uz=1.-0.5*theta*theta;		//Theta is small by assumption
	newU[0]=genrand_pmOne()*uz;
	//Homogenous on the cone
	phi=genrand_real()*2.0*pi;
	newU[1]=sqrt(1.0-uz*uz)*cos(phi);
	newU[2]=sqrt(1.0-uz*uz)*sin(phi);

	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Unrotated orientation (gauss): " );
			pvec( newU,_3D );
		}
	#endif

	//Rotate to be oriented along DIR
	rodriguesRotation( newU,rotAx,rotAngle );
	for( d=0; d<DIM; d++ ) U[d]=newU[d];

	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Rotated orientation (gauss): " );
			pvec( newU,_3D );
		}
	#endif
}

///
/// @brief Generate a random, normalized vector by Maier-Saupe distribution in the large effM*S/KBT limit and 2D.
///
/// This function is currently not in use.
/// It generated a random axis, and performs the Rodrigues rotation algorithm on rotAx using rotAngle. 
/// The result was stored in U[]. 
///
/// @param rotAx The axis around which the rotation will occur. Obtained by cross product of the director field and the x-axis.
/// @param rotAngle The angle of the rotation.
/// @param U The vector components generated with this function.
/// @param KBT The temperature.
/// @param S  The nematic magnitude of the local cell.
/// @param effM The molecular field potential (scaled with dimension).
///
void genrand_maierSaupeGAUSS_2D( double rotAx[],double rotAngle,double U[],double KBT,double S,double effM ) {
	double std,theta;
	double newU[_3D]={0.0};
	int d;

	for( d=0; d<DIM; d++ ) newU[d]=U[d];
	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Old orientation: " );
			pvec( U,_3D );
		}
	#endif

	//Generate a new angle
	std=1./sqrt(2.*effM*S/KBT);
	theta=genrand_gaussGen(0.0,std);
	newU[0]=genrand_pmOne()*cos(theta);
	newU[1]=sin(theta);

	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Unrotated orientation (gauss): " );
			pvec( newU,_3D );
		}
	#endif

	//Rotate to be oriented along DIR
	rodriguesRotation( newU,rotAx,rotAngle );
	for( d=0; d<DIM; d++ ) U[d]=newU[d];

	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Rotated orientation (gauss): " );
			pvec( newU,_3D );
		}
	#endif
}

///
/// @brief Generates a random, normalized vector by Maier-Saupe distribution in the small effM*S/KBT limit and 3D.
///
/// This function uses an expansion around effM*S/KBT to generate a new vector component U.
///
/// @param rotAx The axis around which the rotation will occur. Obtained by cross product of the director field and the x-axis.
/// @param rotAngle The angle of the rotation.
/// @param U The vector components generated with this function.
/// @param KBT The temperature.
/// @param S  The nematic magnitude of the local cell.
/// @param effM The molecular field potential (scaled with dimension)
///
void genrand_maierSaupeEXP_3D( double rotAx[],double rotAngle,double U[],double KBT,double S,double effM ) {
	double BUS,c,R,num,denom,phi,uz;
	double newU[_3D]={0.0};
	int d;

	for( d=0; d<DIM; d++ ) newU[d]=U[d];
	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Old orientation (exp): " );
			pvec( U,_3D );
		}
	#endif

	BUS=effM*S/KBT;
	// The code crashes if a==0 so set a small minimum (if a is small then the fluid is isotropic and shouldn't be run anyway)
	if(BUS<TOL) BUS=TOL;
	//Generate a random number on [0,1+a/3]
	R = genrand_real()*(1.+BUS/3.);
	//Transform R into a random uz
	c=sqrt( 9.*BUS*BUS*BUS*BUS*R*R+4.*BUS*BUS*BUS );
	num=smrtPow( 3.*BUS*BUS*R+c, 1./3. );
	denom=smrtPow( 2.,1./3. );
	uz=num/(denom*BUS) - denom/num;

	newU[0]=genrand_pmOne()*uz;
	//Homogenous on the cone
	phi=genrand_real()*2.0*pi;
	newU[1]=sqrt(1.0-uz*uz)*cos(phi);
	newU[2]=sqrt(1.0-uz*uz)*sin(phi);

	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "BUS=%lf\tnum=%lf\tdenom=%lf\n",BUS,num,denom );
			printf( "phi=%lf\tuz=%lf\n",phi,uz );
			printf( "Unrotated orientation (exp): " );
			pvec( newU,_3D );
		}
	#endif

	//Rotate to be oriented along DIR
	rodriguesRotation( newU,rotAx,rotAngle );
	for( d=0; d<DIM; d++ ) U[d]=newU[d];

	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Rotated orientation (exp): " );
			pvec( newU,_3D );
		}
	#endif
}

///
/// @brief Generate a random, normalized vector by Maier-Saupe distribution using the Metropolis algohritm in 3D.
///
/// This function generates a new vector component U, based on the Metropolis algorithm.
/// The function generates random numbers and see if they are accepted based on the original energy and trial energy (from the random numbers).
///
/// @param DIR The director of the local cell.
/// @param rotAx The axis around which the rotation will occur. Obtained by cross product of the director field and the x-axis.
/// @param rotAngle The angle of the rotation.
/// @param U The vector components generated with this function.
/// @param KBT The temperature.
/// @param S  The nematic magnitude of the local cell.
/// @param effM The molecular field potential (scaled with dimension)
///
void genrand_maierSaupeMetropolis_3D( double DIR[],double rotAx[],double rotAngle,double U[],double KBT,double S,double effM ) {
	double phi,R;
	double newU[_3D]={0.0};
	double uz0,w0,uz1,w1,dw;			//orientation along DIR
	int accept=0;
	int i,d;
	int annealNum;
	int cnt=0;

	for( d=0; d<DIM; d++ ) newU[d]=U[d];
	//Probability only depends on orienation compared to DIR "x"
	uz0 = dotprod( DIR,U,DIM );
	//Metropolis Algorithm
	annealNum=(int)(MCINT+MCSLOPE*effM*S/KBT);
	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Old orientation (metro): " );
			pvec( U,_3D );
			printf( "uz0= %lf\n",uz0 );
			printf( "annealNum= %d\n",annealNum );
		}
	#endif
	for( i=0; i<annealNum; i++ ) {
		accept=0;
		//Calculate original energy and trial energy --- don't include constant (wrt u.n) term
		w0 = effM*S*uz0*uz0/KBT;
		uz1 = genrand_real();			//Try a new uz
		w1 = effM*S*uz1*uz1/KBT;
		dw=w1-w0;
		//If the new energy is less than the old energy then accept
		if( dw>=0. ) accept=1;
		//Else generate a new random number
		else {
			R = genrand_real();		//Random number for "uphill move"
			if( R<=exp(dw) ) accept=1;
		}
		if(accept) {
			//We have already generated the random vector relative to the "z-axis"
			cnt+=1;
			uz0=uz1;
		}
	}
	//Generate the other dimensions based on the updated uz
	//Homogenous on the cone
	phi=genrand_real()*2.0*pi;
	newU[0]=genrand_pmOne()*uz0;
	newU[1]=sqrt(1.0-uz0*uz0)*cos(phi);
	newU[2]=sqrt(1.0-uz0*uz0)*sin(phi);
	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Unrotated orientation (metro): " );
			pvec( newU,_3D );
			printf( "Number of accepted moves (metro): %d\n",cnt );
			//NOTICE: When fails it is ALWAYS zero accepted moves!!!!
		}
	#endif

	if( cnt>0 ) {
		rodriguesRotation( newU,rotAx,rotAngle );
		for( d=0; d<DIM; d++ ) U[d]=newU[d];
	}
	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Rotated orientation (metro): " );
			pvec( newU,_3D );
		}
	#endif
}

///
/// @brief Generate a random, normalized vector by Maier-Saupe distribution using the Metropolis algohritm in 2D.
///
/// This function generates a new vector component U, based on the Metropolis algorithm.
/// The function generates random numbers and see if they are accepted based on the original energy and trial energy (from the random numbers).
/// In 2D must consider the distribution of the ANGLE between u and DIR
///
/// @param DIR The director of the local cell.
/// @param rotAx The axis around which the rotation will occur. Obtained by cross product of the director field and the x-axis.
/// @param rotAngle The angle of the rotation.
/// @param U The vector components generated with this function.
/// @param KBT The temperature.
/// @param S  The nematic magnitude of the local cell.
/// @param effM The molecular field potential (scaled with dimension).
///
void genrand_maierSaupeMetropolis_2D( double DIR[],double rotAx[],double rotAngle,double U[],double KBT,double S,double effM ) {
	double newU[_3D]={0.0};
	double R,ang0,uz0,w0,ang1,uz1,w1,dw;			//orientation along DIR
	int i,d,annealNum;
	int accept=0;
	int cnt=0;

	for( d=0; d<DIM; d++ ) newU[d]=U[d];
	//Probability only debends on orienation compared to DIR --- angle in 2D
	ang0 = signedAngle( DIR,U,DIM );
	//NOTICE: If a mesogen initializes with an angle pi/2 then it will always stay at that initial angle and there will be an artificial spike in the orientation angle distribution.
	if( feq(ang0,0.5*pi) || feq(ang0,-0.5*pi) ) ang0 = pi*(0.5-genrand_real());
	uz0 = dotprod( DIR,U,DIM );
	//Metropolis Algorithm
	annealNum=(int)(MCINT+MCSLOPE*effM*S/KBT);
	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Old orientation: " );
			pvec( U,_3D );
			printf( "ang0= %lf\n",ang0 );
			printf( "annealNum= %d\n",annealNum );
		}
	#endif
	for( i=0; i<annealNum; i++ ) {
		accept=0;
		//Calculate original energy and trial energy --- don't include constant (wrt u.n) term
		w0 = effM*S*uz0*uz0/KBT;
		ang1 = pi*(0.5-genrand_real());		//Try a new angle
		uz1 = cos(ang1);
		w1 = effM*S*uz1*uz1/KBT;
		dw=w1-w0;
		//If the new energy is less than the old energy then accept
		if( dw>=0. ) accept=1;
		//Else generate a new random number
		else {
			R = genrand_real();		//Random number for "uphill move"
			if( R<=exp(dw) ) accept=1;
		}
		if(accept) {
			//We have already generated the random vector relative to the "z-axis"
			cnt+=1;
			//Update
			ang0=ang1;
			uz0=uz1;
		}
	}
	//Save the updated angle
	newU[0]=uz0;
	newU[1]=sin(ang0);

	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Unrotated orientation (metro): " );
			pvec( newU,_3D );
			printf( "Number of accepted moves (metro): %d\n",cnt );
			//NOTICE: When fails it is ALWAYS zero accepted moves!!!!
		}
	#endif

	if( cnt>0 ) {
		rodriguesRotation( newU,rotAx,rotAngle );
		for( d=0; d<DIM; d++ ) U[d]=newU[d];
	}

	#ifdef DBG
		if( DBUG == DBGLCCOL ) {
			printf( "Rotated orientation (metro): " );
			pvec( newU,_3D );
		}
	#endif
}

///
/// @brief Perform the multi-particle orientation collision operation on a single cell.
///
/// Does the multi-particle orientation collision in a single cell.
/// Collision has randomly sampled orientations (which obey the Maier-Saupe distribution) around the local director.
///
/// @param CL Class containing cell data (pointed to local cell considered).
/// @param SP List of species.
/// @param KBT The temperature.
/// @param MFPOT The molecular field potential.
/// @param dt The timestep.
/// @param SG The global nematic order value.
/// @param LC Variable to define if Nematic LC using an isotropic fluid (0), the local nematic order  value (1) or global nematic order value (2).
///
void LCcollision( cell *CL,spec *SP,double KBT,double MFPOT,double dt,double SG,int LC ) {
	int i,id;
	double DIR[_3D]={0.0},dU[_3D]={0.0};	//The director, the difference in orientation
	double rotAx[_3D],xaxis[_3D]={0.0},rotAngle;
	double rfc;			//Rotational friction coefficient
	double S;			//Local scalar order parameter
	double MFPOT_scaled=MFPOT*0.5*(double)DIM;
	particleMPC *tmpc;		//Temporary particleMPC

	if( LC==LCL ) S=CL->S;		//Use the local scalar order parameter in collision
	else if( LC==LCG ) S=SG;	//Use the global scalar order parameter in collision
	else{
		printf( "Error: Unexpected value of LC=%d in LCcolloision().\n",LC );
		exit( 1 );
	}

	// Set arrays
	xaxis[0]=1.0;
	for( i=0; i<DIM; i++ ) DIR[i] = CL->DIR[i];
	#ifdef DBG
		if( DBUG==DBGLCCOL || DBUG==DBGESCAPE ) {
			printf( "Director: " );
			pvec( DIR,_3D );
		}
	#endif
	#ifdef DBG
		if( DBUG>=DBGLCCOL ) if( dotprod(DIR,DIR,_3D)<=TOL ) {
			printf( "WARNING: Director is zero n=" );
			pvec( DIR,_3D );
			printf( "\tS=%lf\n",S );
			printf( "\tCell pop=%i\n",CL->POPSRD );
			tmpc = CL->pp;
			while( tmpc!=NULL ) {
				id = tmpc->SPID;
				pcoord( *tmpc );
				//Increment link in list
				tmpc = tmpc->next;
			}
		}
	#endif
	//Find the axis that the x-axis must be rotated about
	crossprod( DIR,xaxis,rotAx );
	rotAngle = -1.*absAngle( DIR,xaxis,DIM );	//The minus sign seems to be needed in 2D. I haven't check 3D yet!!!
	norm( rotAx,_3D );
	#ifdef DBG
		if( DBUG==DBGLCCOL || DBUG==DBGESCAPE ) {
			printf( "Rotation axis: " );
			pvec( rotAx,_3D );
			printf( "Rotation angle: %lf\n",rotAngle );
		}
	#endif

	//Generate random orientations for MPC particles
	if( CL->POPSRD>1 ) {
		//Collision of MPC particles
		tmpc = CL->pp;
		while( tmpc!=NULL ) {
			id = tmpc->SPID;
			rfc=(SP+id)->RFC;
			//Zero torque
			for( i=0; i<_3D; i++ ) tmpc->T[i] = 0.;
			//Save old orientation
			for( i=0; i<DIM; i++ ) dU[i]=tmpc->U[i];
			//Metropolis generation of new orientation and record the change
			#ifdef DBG
				if( DBUG==DBGLCCOL || DBUG==DBGESCAPE ) {
					printf( "---\nMaier-Suape\n" );
					printf( "MPC particle orientation:" );
					pvec( tmpc->U,DIM );
				}
			#endif
			genrand_maierSaupe( DIR,rotAx,rotAngle,tmpc->U,KBT,S,MFPOT_scaled );
			#ifdef DBG
				if( DBUG==DBGLCCOL || DBUG==DBGESCAPE ) {
					printf( "New orientation: " );
					pvec( tmpc->U,DIM );
				}
			#endif
			//Calculate rate of change of orientation
			for( i=0; i<DIM; i++ ) dU[i]=(tmpc->U[i]-dU[i])/dt;
			//Calculate angular velocity (save in torque for now)
			//--- must have 3rd dimension even if 2D
			crossprod( tmpc->U,dU,tmpc->T );
			//Turn angular velocity into torque
			for( i=0; i<_3D; i++ ) tmpc->T[i] *= rfc;
			//Increment link in list
			tmpc = tmpc->next;
		}
		#ifdef DBG
			if( DBUG==DBGLCCOL || DBUG==DBGESCAPE ) printf( "\n" );
		#endif
	}
}

///
/// @brief Apply a torque on the cell director based on velocity gradients.
///
/// The velocity gradient shears the rods in the cell, applying a torque on them which affects the director of the cell. It is assumed there is no external torque so
/// the rod rotates due to the shear (hydrodynamical torque).
/// See Jeffery, G.B.: The motion of ellipsoidal particles immersed in a viscous flow. Proc. R. Soc. A 102, 161 (1922), https://doi.org/10.1098/rspa.1922.0078.
/// Currently uses the Larson rotational rate (see 'The structure and rheology of complex fluids', Larson, p. 280, eq. 6-26).
///
/// @param CL Class containing cell data (pointed to local cell considered).
/// @param SP List of species.
/// @param dt The integrator timestep.
/// @see larsonRotRate()
///
void jefferysTorque( cell *CL,spec *SP,double dt ) {
	int i,id;
	double CHIHI;			// Susceptibility to shear
	double tumble;			//Tumbling parameter of nematogens
	double dudt[_3D],u[_3D],w[_3D];
	particleMPC *tmpc;		//Temporary particleMPC
	double theta;

	//Apply the shearing torque
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		tumble=(SP+id)->TUMBLE;
		CHIHI=(SP+id)->CHIHI;
		//Determine the rate of change of orientation
		for( i=0; i<_3D; i++ ) u[i]=tmpc->U[i];
		#ifdef DBG
			if( DBUG == DBGJEFF ) {
				printf("u\t");
				pvec( u,_3D );
				printf("|u|\t%lf\n",length(u,_3D));
			  }
		#endif
		//Calculate the rotation rate
		larsonRotRate(dudt,w,u,CL->E,tumble);
		//brielsRotRate(dudt,w,u,CL->E);
		//saintillanRotRate(dudt,w,u,CL->E,tumble);
		//Scale with CHIHI ---"Shear Susceptibility" is heuristically less than 1 and theoretically CHIHI=1.
		for( i=0; i<_3D; i++ ) w[i] *= CHIHI;
		//Calculate new orientation from the angular velocity
		for( i=0; i<_3D; i++ ) w[i]*=dt;
		theta=length(w,_3D);
		norm(w,_3D);
		//Do the rotation
		rodriguesRotation( tmpc->U,w,theta );
		#ifdef DBG
			if( DBUG == DBGJEFF ) {
				printf("theta=%lf\n",theta*180./pi);
				printf("U\t");
				pvec( tmpc->U,_3D );
				printf("|U|\t%lf\n",length(tmpc->U,_3D));
				printf("\n");
			  }
		#endif

		//Increment link in list
		tmpc = tmpc->next;
	}
}

///
/// @brief Calculate the rotation rate of the rod (dudt) and angular velocity (w) using Larson.
///
/// 'The structure and rheology of complex fluids', Larson, p. 448, eq. 10-3 used for the Jefferey torque.
/// Would be Larson p. 280, eq. 6-26 if the tumbling parameter was tumbleParam=(P*P-1.)/(P*P+1.) where P=aspect ratio.
/// Velocity gradient tensor E --- first index [i] is on velocity, second [j] on derivative --- E[i][j]= dv[i]/dx[j].
///
/// @param dudt The rotation rate.
/// @param w The angular velocity.
/// @param u The orientation of the rod.
/// @param E The velocity gradient.
/// @param tumbleParam The tumbling parameter.
///
void larsonRotRate(double dudt[],double w[],double u[],double E[_3D][_3D],double tumbleParam) {
//Loops done explicitly by hand to save computational time

/***FIXME: the below code, while meant to be faster, doesn't given the same results as the original code.***/

// 	int i,j;
// 	double T1,T2,t2const,D[_3D][_3D],W[_3D][_3D];

// 	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) {
// 			D[j][i] = (E[i][j]+E[j][i]);
// 			W[j][i] = (E[i][j]-E[j][i]);
// 	}
// 	//Calculate rate of rotation of orientation
// 	//Let j be the vector index
//     if(DIM==_3D) {
//         //This term is the same for every j. So calculate outside of loop
//         t2const = u[0]*(u[0]*D[0][0]+u[1]*D[1][0]+u[2]*D[2][0]) + u[1]*(u[0]*D[0][1]+u[1]*D[1][1]+u[2]*D[2][1]) + u[2]*(u[0]*D[0][2]+u[1]*D[1][2]+u[2]*D[2][2]);
//         for( j=0; j<_3D; j++ ) {
//             //Vorticity term --- u dotted into rate of antisymmetric deformation tensor
//             T1 = u[0]*W[0][j]+u[1]*W[1][j]+u[2]*W[2][j];
//             //"Second" strain rate term --- Third order uuu double dotted into rate of deformation tensor
//             T2 = -u[j]*t2const;
//             //"First" strain rate term --- u dotted into rate of symmetric deformation tensor
//             T2 += u[0]*D[0][j]+u[1]*D[1][j]+u[2]*D[2][j];
//             //The factor of half comes from not including them in every vorticity/shear rate
//             dudt[j]=0.5*(T1+tumbleParam*T2);
//         }
// 	}
// 	else if(DIM==_2D) {
//         //This term is the same for every j. So calculate outside of loop
//         t2const=u[0]*(u[0]*D[0][0]+u[1]*D[1][0]) +  u[1]*(u[0]*D[0][1]+u[1]*D[1][1]);
//         //Same as 3D but explicitly removes z-components
//         for( j=0; j<_2D; j++ ) {
//             //Vorticity term --- u dotted into rate of antisymmetric deformation tensor
//             T1 = u[0]*W[0][j]+u[1]*W[1][j];
//             //"Second" strain rate term --- Third order uuu double dotted into rate of deformation tensor
//             T2 = -u[j]*t2const;
//             //"First" strain rate term --- u dotted into rate of symmetric deformation tensor 44
//             T2 += u[0]*D[0][j]+u[1]*D[1][j];
//             //The factor of half comes from not including them in every vorticity/shear rate
//             dudt[j]=0.5*(T1+tumbleParam*T2);
//         }
// 		//Set the last component to zero. Might not be necessary but dudt later goes into crossprod() with 3D vectors
// 		dudt[2]=0.0;
// 	}
// 	else printf("Warning: larsonRotRate() only programmed for DIM={3,2}, not DIM=%d\n",DIM);

// 	//Convert to angular velocity
// 	crossprod( u,dudt,w );
// 	#ifdef DBG
// 		if( DBUG == DBGJEFF ) {
// 			printf("dudt\t");
// 			pvec( dudt,_3D );
// 		  }
// 	#endif
// }
// void larsonRotRateOLD_AND_SLOW(double dudt[],double w[],double u[],double E[_3D][_3D],double tumbleParam) {
int i,j,k;
	double uuuD[_3D],uw[_3D],uD[_3D];
	double D[_3D][_3D],W[_3D][_3D];

	for( i=0; i<_3D; i++ ) {
		uuuD[i]=0.0;
		uw[i]=0.0;
		uD[i]=0.0;
		for( j=0; j<_3D; j++ ) {
			D[j][i] = 0.5*( E[i][j] + E[j][i] );
			W[j][i] = 0.5*( E[i][j] - E[j][i] );
		}
	}
	//Calculate rate of rotation of orientation
	//Let j be the vector index
	for( j=0; j<_3D; j++ ) {
		for (i = 0; i < _3D; i++){
			//Third order uuu double dotted into rate of deformation tensor = vector
			for( k=0; k<_3D; k++ ) uuuD[j] += u[j]*( u[i]*D[i][k]*u[k] );
			//Vector u dotted into rate of symmetric deformation tensor = vector
			uD[j] += u[i]*D[i][j];
			//Vector u dotted into rate of antisymmetric deformation tensor = vector
			uw[j] += u[i]*W[i][j];
		}
		dudt[j]=uw[j]+tumbleParam*(uD[j]-uuuD[j]);
	}

	//Convert to angular velocity
	crossprod( u,dudt,w );
	#ifdef DBG
		if( DBUG == DBGJEFF ) {
			printf("uw\t");
			pvec( uw,_3D );
			printf("l*(uD-uuuD)\t");
			printf("(%lf,%lf,%lf)\n",dudt[0]-uw[0],dudt[1]-uw[1],dudt[2]-uw[2]);
			printf("dudt\t");
			pvec( dudt,_3D );
		  }
	#endif
}

///
/// @brief Calculate the rotation rate of the rod (dudt) and angular velocity (w) using Dhont and Briels.
///
/// Currently not in use.
/// Calculates the rotation rate of the rod (dudt) and angular velocity (w) according to Dhont and Briels,
/// 'Rod-like Brownian Particles in Shear Flow' ( https://onlinelibrary.wiley.com/doi/10.1002/9783527617067.ch3a ), for the Jefferey torque.
/// Velocity gradient tensor E --- first index [i] is on velocity, second [j] on derivative --- E[i][j]= dv[i]/dx[j].
///
/// @param dudt The rotation rate.
/// @param w The angular velocity.
/// @param u The rod orientation.
/// @param E The velocity gradient.
///
void brielsRotRate(double dudt[],double w[],double u[],double E[_3D][_3D]) {
	double t1[_3D];			//Temporary vectors that take various values --- Pay attention
	//Calculate rate of rotation of orientation
	dotprodMatVec( E,u,t1,_3D );	//t1=dot( E,u )
	crossprod( u,t1,w );	//w=cross( u,dot( E,u ) )
	//Convert to rotation rate
	crossprod( w,u,dudt );	//dudt=cross( w,u )
	#ifdef DBG
		if( DBUG == DBGJEFF ) {
			printf("E dot u\t");
			pvec( t1,_3D );
			printf("w\t");
			pvec( w,_3D );
			printf("dudt\t");
			pvec( dudt,_3D );
		  }
	#endif
}

///
/// @brief Calculate the rotation rate of the rod (dudt) and angular velocity (w) using Saintillan.
///
/// Currently not in use.
/// Calculates the rotation rate of the rod (dudt) and angular velocity (w) according to Saintillan review ( https://doi.org/10.1146/annurev-fluid-010816-060049 ).
/// Velocity gradient tensor E --- first index [i] is on velocity, second [j] on derivative --- E[i][j]= dv[i]/dx[j].
///
/// @param dudt The rotation rate.
/// @param w The angular velocity.
/// @param u The rod orientation.
/// @param E The velocity gradient.
/// @param tumbleParam The tumbling parameter.
///
void saintillanRotRate(double dudt[],double w[],double u[],double E[_3D][_3D],double tumbleParam) {
	int i,j;
	double t1[_3D],t2[_3D],Ipp[_3D][_3D];

	//Calculate rate of rotation of orientation
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) Ipp[i][j] = -u[i]*u[j];
	for( i=0; i<_3D; i++ ) Ipp[i][i]+=1.;

	// split E into symmetric and antisymmetric parts
	double D[_3D][_3D],W[_3D][_3D];
	for( i=0; i<_3D; i++ ) {
		for( j=0; j<_3D; j++ ) {
			D[j][i] = 0.5*( E[i][j] + E[j][i] );
			W[j][i] = 0.5*( E[i][j] - E[j][i] );
		}
	}
	dotprodMatVec( D,u,t1,_3D );
	dotprodMatVec( W,u,t2,_3D );
	for (i=0; i<_3D; i++) t1[i] = tumbleParam*t1[i] + t2[i];

	dotprodMatVec( Ipp,t1,dudt,_3D );
	//Convert to angular velocity
	crossprod( u,dudt,w );
	#ifdef DBG
		if( DBUG == DBGJEFF ) {
			printf("dudt\t");
			pvec( dudt,_3D );
			printf("w\t");
			pvec( w,_3D );
		  }
	#endif
}

///
/// @brief Calculate the torque and apply rotation due to the external magnetic field on an MPCD particle.
///
/// Includes a check to see if the angle "swings past" the magnetic field (shouldn't be allowed by overdamped assumption).
///
/// @param pMPC Pointer to specific MPCD particle.
/// @param SP List of species.
/// @param dt The rod orientation.
/// @param MAG The external applied field.
///
void magTorque( particleMPC *pMPC,spec *SP,double dt,double MAG[] ) {
	int i,id;
	double uH,magMAG;
	double chia,rfc,theta,theta0;
	double mT[_3D],U[_3D],angVel[_3D],rotAxis[_3D];

	id = pMPC->SPID;
	chia = (SP+id)->CHIA;
	rfc = (SP+id)->RFC;
	for( i=0; i<_3D; i++ ) 	mT[i]=0.;
	for( i=0; i<DIM; i++ ) U[i]=pMPC->U[i];
	for( i=DIM; i<_3D; i++ ) U[i] = 0.;
	magMAG = length( MAG,DIM );

	//Calculate torque:
	//Calculate dot product of U and magnetic field
	uH=dotprod( U,MAG,DIM );
	//Calculate the original angle between U and magnetic field
	theta0=acos( uH/magMAG );
	//Calculate the cross product of U and magnetic field (which is in the direction of the torque)
	crossprod( U,MAG,mT );
	#ifdef DBG
		if( DBUG == DBGMAG ) {
			printf( "Test cross product:\n" );
			printf( "Old:\t" );
			crossprod( U,MAG,mT );
			pvec( mT,_3D );
			printf( "New:\t" );
			crossprod( U,MAG,mT );
			pvec( mT,_3D );
		}
	#endif
	//Calculate the torque
	for( i=0; i<_3D; i++ ) mT[i] *= chia*uH;
	//Update
	//Add this torque to the total torque
	for( i=0; i<_3D; i++ ) pMPC->T[i]+=mT[i];	//Torque must always be 3D
	//Calculate change in orientation due to magnetic torque
	//Angular velocity
	for( i=0; i<_3D; i++ ) angVel[i]=mT[i]/rfc;
	//Angular velocity axis (just normalized)
	for( i=0; i<_3D; i++ ) rotAxis[i]=angVel[i];
	norm( rotAxis,_3D );
	//Angular speed times time gives the change in angle about the rotation axis
	theta=length(angVel,_3D)*dt;
	if( theta>=theta0 ) for( i=0; i<DIM; i++ ) pMPC->U[i] = MAG[i]/magMAG;
	else {
		//Do the rotation
		rodriguesRotation( pMPC->U,rotAxis,theta );
	}
	#ifdef DBG
		if( DBUG == DBGMAG ) {
			printf("aChi = %lf\t friction = %lf\n",chia,rfc);
			printf("MAG\t");
			pvec( MAG,_3D );
			printf("U\t");
			pvec( U,_3D );
			printf("u dot H = %lf\n",uH);
			printf("theta0=%lf\n",theta0*180./pi);
			printf("Torque\t");
			pvec( mT,_3D );
			printf("angular velocity\t");
			pvec( angVel,_3D );
			printf("rotation axis\t");
			pvec( rotAxis,_3D );
			printf("theta=%lf\n",theta*180./pi);
			printf("U'\t");
			pvec( pMPC->U,_3D );
			printf("\n");
		  }
	#endif
}

///
/// @brief Apply the torque due to external magnetic field on all MPCD particles.
///
/// @param pp pointer to first MPCD particle.
/// @param SP List of species.
/// @param dt The rod orientation.
/// @param MAG The external applied field.
///
void magTorque_all( particleMPC *pp,spec *SP,double dt,double MAG[] ) {
	int i;
	for( i=0; i<GPOP; i++ ) magTorque( (pp+i),SP,dt,MAG );
}

///
/// @brief  Rotate all the orientations in a given cell towards the external magnetic field.
///
/// Currently not in use.
///
/// @param CL Class containing cell data (pointed to local cell considered).
/// @param SP List of species.
/// @param dt The rod orientation.
/// @param MAG The external applied field.
///
void magTorque_CL( cell *CL,spec *SP,double dt,double MAG[] ) {
	int i,id;
	double nH;
	double chia,rfc,theta;
	double mT[_3D],DIR[_3D],w[_3D];
	particleMPC *tmpc;		//Temporary particleMPC

	// Set arrays
	for( i=0; i<DIM; i++ ) DIR[i] = CL->DIR[i];
	for( i=DIM; i<_3D; i++ ) DIR[i] = 0.;
	for( i=0; i<_3D; i++ ) 	mT[i]=0.;
	#ifdef DBG
		if( DBUG == DBGMAG ) {
			printf( "Director: " );
			pvec( DIR,_3D );
		}
	#endif
	//Calculate torque for the entire cell
	nH=dotprod( DIR,MAG,DIM );
	crossprod( DIR,MAG,mT );
	for( i=0; i<_3D; i++ ) mT[i] *= nH;		//Will need to multiply chia

	//Rotate each MPC particles due to the torque on the whole cell
	if( CL->POPSRD>1 ) {
		//Collision of MPC particles
		tmpc = CL->pp;
		while( tmpc!=NULL ) {
			id = tmpc->SPID;
			chia = (SP+id)->CHIA;
			rfc = (SP+id)->RFC;
			//Update torque on MPC particles
			for( i=0; i<_3D; i++ ) tmpc->T[i]+=chia*mT[i];	//Torque must always be 3D (needs to be multiplied by particle's susceptibility)
			//Calculate change in orientation due to magnetic torque
			//Angular velocity
			for( i=0; i<_3D; i++ ) w[i]=chia*mT[i]/rfc;
			//Angular velocity times time gives the angle change
			for( i=0; i<_3D; i++ ) w[i]*=dt;
			theta=length(w,_3D);
			norm(w,_3D);
			//Do the rotation
			rodriguesRotation( tmpc->U,w,theta );
			#ifdef DBG
				if( DBUG == DBGMAG ) {
					printf("theta=%lf\n",theta*180./pi);
					printf("DIR\t");
					pvec( DIR,_3D );
					printf("|DIR|\t%lf\n",length(DIR,_3D));
					printf("w\t");
					pvec( w,_3D );
					printf("U'\t");
					pvec( tmpc->U,_3D );
					printf("|U'|\t%lf\n",length(tmpc->U,_3D));
					printf("\n");
				  }
			#endif
			//Increment link in list
			tmpc = tmpc->next;
		}
		#ifdef DBG
			if( DBUG == DBGMAG ) printf( "\n" );
		#endif
	}
}

///
/// @brief Find the average global scalar nematic order of all MPCD particles.
///
/// Constructs a tensor S, which is filled with the calculated values of the nematic orientation tensor Q. This S tensor is then overwritten with its
/// own eigenvectors, and the largest eigenvalue is returned as the order parameter; in 3-D, a correction is added with the following two eigenvalues.
/// If the largest eigenvalue is negative, a warning is issued
/// and the corresponding eigenvector is printed.
///
/// @param p Pointer to first MPCD particle.
/// @param LC Variable to define if Nematic LC using an isotropic fluid (0), the local nematic order  value (1) or global nematic order value (2).
/// @param avDIR The average director of all the MPCD particles in the simulation calculated with this routine.
/// @return The global scalar order parameter of the MPCD particles.
///
double avOrderParam( particleMPC *p,int LC,double avDIR[] ) {
	int i,j,k;
	double **S,eigval[_3D],avS;
	double U[_3D];
	double GPOPinv,fDIM,c;

	avS=0.;
	GPOPinv=1./((double)GPOP);
	fDIM=(double)DIM;
	c=1./(fDIM-1.);

	//Calculate the tensor order parameter
	//Allocate memory for S and zero
	S = malloc ( DIM * sizeof( *S ) );
	for( i=0; i<DIM; i++ ) S[i] = malloc ( DIM * sizeof( *S[i] ) );
	for( j=0; j<DIM; j++ ) for( k=0; k<DIM; k++ ) S[j][k] = 0.0;
	for( j=0; j<_3D; j++ ) U[j] = 0.0;

	//This uses fact that S is always symmetric by only calculating half then copying
	if( LC!=ISOF ) for( i=0; i<GPOP; i++ ) {
		for( j=0; j<DIM; j++ ) U[j] = (p+i)->U[j];
		for( j=0; j<DIM; j++ ) for( k=j; k<DIM; k++ ) S[j][k] += U[j] * U[k];
	}
	else for( i=0; i<GPOP; i++ ) {
		normCopy( (p+i)->V,U,DIM );
		for( j=0; j<DIM; j++ ) for( k=j; k<DIM; k++ ) S[j][k] += U[j] * U[k];
	}
	// for( j=0; j<DIM; j++ ) for( k=j; k<DIM; k++ ) S[j][k] *= GPOPinv;
	// for( i=0; i<DIM; i++ ) S[i][i] -= 1/fDIM;
	// for( j=0; j<DIM; j++ ) for( k=j; k<DIM; k++ ) S[j][k] *= c*fDIM;
	// //Fill in other half
	// for( j=0; j<DIM; j++ ) for( k=j+1; k<DIM; k++ ) S[k][j] = S[j][k];
	// Make sum into average
	for( i=0; i<DIM; i++ ) for( j=i; j<DIM; j++ ) S[i][j] *= c*fDIM*GPOPinv;
	for( i=0; i<DIM; i++ ) S[i][i] -= c;
	// Copy other half
	for( i=0; i<DIM; i++ ) for( j=i+1; j<DIM; j++ ) S[j][i] = S[i][j];
	//Solve the eigensystem to find the order parameter as the largest eigenvalue
	solveEigensystem( S,DIM,eigval );
	//The scalar order parameter is the largest eigenvalue
	//The scalar order parameter is the largest eigenvalue which is given first ie eigval[0]
	// But can be better approximated (cuz fluctuates about zero) by -2* either of the negative ones (or the average)
	if(DIM==_3D) avS = -1.*(eigval[1]+eigval[2]);
	else avS=eigval[0];
	if( avS<1./(1.-DIM) ){
		printf("Warning: Global scalar order parameter <0\n");
		printf("Eigenvalues=");
		pvec(eigval,DIM);
		printf("Eigenvectors=");
		for( j=0; j<DIM; j++ ) pvec(S[j],DIM);
	}
	// The director is the eigenvector corresponding to the largest eigenvalue
	for( i=0; i<DIM; i++ ) avDIR[i] = S[0][i];

	for( i=0; i<DIM; i++ ) free( S[i] );
	free( S );
	return avS;
}

///
/// @brief Find and return the average global scalar fourth order mode of all MPCD particles.
///
/// It uses the already calculated average director (DIR) of the MPCD particles. Uses Legendre polynomials in 3D, Chebyshev in 2D.
///
/// @param p Pointer to first MPCD particle.
/// @param LC Variable to define if Nematic LC using an isotropic fluid (0), the local nematic order  value (1) or global nematic order value (2).
/// @param DIR The average director of all the MPCD particles in the simulation.
/// @return The average global scalar fourth order mode.
///
double avS4( particleMPC *p,int LC,double DIR[] ) {
	double U[_3D];
	double un,un2,un4;	//orientation dot director
	double S4;
	int i;

	//Initialize
	un2 = 0.;
	un4 = 0.;
	S4 = 0.;
	//Calculate average un2 and un4
	for( i=0; i<GPOP; i++ ) {
		if( LC!=ISOF ) normCopy( (p+i)->U,U,DIM );
		else normCopy( (p+i)->V,U,DIM );
		un = dotprod( U,DIR,DIM );
		//Sum
		un2 += un*un;
		un4 += un*un*un*un;
	}
	un2 /= ((double) GPOP);
	un4 /= ((double) GPOP);
	//If 3D use Legendre polynomial; if 2D use Chebychev
	if( DIM==_3D ) S4 = 0.125*(35.*un4-30.*un2+3.);
	else if( DIM==_2D ) S4 = 8.*un4-8.*un2+1.;
	return S4;
}

///
/// @brief Calculate the tensor order parameter from a set of orientation vectors or MPCD particles.
///
/// @param pMPC A pointer to the first MPCD particle considered.
/// @param S The calculated tensor order parameter (taken as the sum, needs to be normalized).
///
void addToTensOrderParam( particleMPC *pMPC,double **S ) {
	int i,j;

	while(pMPC!=NULL) {
		//This uses fact that S is always symmetric by only calculating half then copying
		for( i=0; i<DIM; i++ ) for( j=i; j<DIM; j++ ) S[i][j] += pMPC->U[i] * pMPC->U[j];
		//Increment link in list
		pMPC = pMPC->next;
	}
}

///
/// @brief Calculate the sum of the tensor order parameter from a set of <b>velocity</b> vectors of MPCD particles.
///
/// Used for isotropic fluids where velocity is same as director field.
///
/// @param pMPC A pointer to the first MPCD particle considered.
/// @param S The calculated tensor order parameter (taken as the sum, needs to be normalized).
///
void addToTensOrderParamVel( particleMPC *pMPC,double **S ) {
	int i,j;
	double U[_3D];

	while(pMPC!=NULL) {
		//This uses fact that S is always symmetric by only calculating half then copying
		normCopy( pMPC->V,U,DIM );
		for( i=0; i<DIM; i++ ) for( j=i; j<DIM; j++ ) S[i][j] += U[i] * U[j];
		//Increment link in list
		pMPC = pMPC->next;
	}
}

///
/// @brief Calculate the tensor order parameter in a cell.
///
/// For isotropic fluids (LC=0) velocity is same as director field, otherwise normal director field.
///
/// @param CL Class containing cell data (pointed to local cell considered).
/// @param S The calculated tensor order parameter.
/// @param LC Variable to define if Nematic LC using an isotropic fluid (0), the local nematic order  value (1) or global nematic order value (2).
///
void tensOrderParam( cell *CL,double **S,int LC ) {
	int i,j;
	particleMPC *pMPC;	//Temporary pointer to MPC particles
	double POPinv,fDIM,c;

	fDIM = (double) DIM;
	c=1./(fDIM-1.);
	//Zero the order parameter
	for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) S[i][j] = 0.0;
	//Calculate the order parameter tensor
	if( CL->pp!=NULL ) {
		pMPC = CL->pp;
		POPinv = 1./((double) CL->POPSRD);
		//Calculate the order parameter tensor in this cell
		if( LC!=ISOF ) addToTensOrderParam( pMPC,S );
		else addToTensOrderParamVel( pMPC,S );

		// Make sum into average
		for( i=0; i<DIM; i++ ) for( j=i; j<DIM; j++ ) S[i][j] *= c*fDIM*POPinv;
		for( i=0; i<DIM; i++ ) S[i][i] -= c;
		// Copy other half
		for( i=0; i<DIM; i++ ) for( j=i+1; j<DIM; j++ ) S[j][i] = S[i][j];
	}
}

///
/// @brief Calculate the tensor order parameter in a cell and its immediate surrounding.
///
/// This routine calculates the tensor order parameter in a cell, its neighbours, and their diagonal next-nearest neighbours.
/// For isotropic fluids (LC=0) velocity is same as director field, otherwise normal director field.
///
/// @param CL Class containing cell data.
/// @param S The calculated tensor order parameter.
/// @param LC Variable to define if Nematic LC using an isotropic fluid (0), the local nematic order  value (1) or global nematic order value (2).
/// @param a x index of cell considered.
/// @param b y index of cell considered.
/// @param c z index of cell considered.
///
void tensOrderParamNNN( cell ***CL,double **S,int LC,int a,int b,int c ) {
	int i,j,x,y,z;
	int POP;
	particleMPC *pMPC;	//Temporary pointer to MPC particles
	double POPinv,fDIM,cnst;

	//Initialize
	fDIM = (double) DIM;
	cnst=1./(fDIM-1.);
	POP=0;
	for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) S[i][j] = 0.0;

	//Add to the order parameter tensor from itself, its neighbours and it's diagonal next-nearest neighbours
	for( x=a-1; x<=a+1; x++ ) for( y=b-1; y<=b+1; y++ ) for( z=c-1; z<=c+1; z++ ) {
		if( ( x>=0 && y>=0 && z>=0 ) && ( x<=XYZ[0] && y<=XYZ[1] && z<=XYZ[2] ) ) {
			if( CL[x][y][z].pp!=NULL ) {
				pMPC = CL[x][y][z].pp;
				POP += CL[x][y][z].POPSRD;
				//Calculate the order parameter tensor in this cell
				if( LC!=ISOF ) addToTensOrderParam( pMPC,S );
				else addToTensOrderParamVel( pMPC,S );
			}
		}
	}
	//Average, etc
	if( POP > 0 ) {
		POPinv = 1./((double) POP);
		// Make sum into average
		for( i=0; i<DIM; i++ ) for( j=i; j<DIM; j++ ) S[i][j] *= cnst*fDIM*POPinv;
		for( i=0; i<DIM; i++ ) S[i][i] -= cnst;
		// Copy other half
		for( i=0; i<DIM; i++ ) for( j=i+1; j<DIM; j++ ) S[j][i] = S[i][j];
	}
}

///
/// @brief Find and return the Binder cumulant for a bin size L of cells.
///
/// The Binder cumulant is calculated over the orientation parameter (velocity if isotropic fluid).
/// See https://en.wikipedia.org/wiki/Binder_parameter.
///
/// @param CL All of the MPCD cells (including the linked list of particles in each cell). 
/// @param L The bin size considered.
/// @param LC Variable to define if Nematic LC using an isotropic fluid (0), the local nematic order  value (1) or global nematic order value (2).
/// @return The Binder cumulant value.
///
double binderCumulant( cell ***CL,int L,int LC ) {
	int nBins[_3D]={0.0};			//Number of cumulant bins
	int i,j,k,a,b,c,d;			//Indices
	double UL,thisS,avS,avS2,avS4;		//Average of the square and power 4 of order parameter and the Binder cumulant
	double **S,eigval[_3D];			//Order parameter tensor and it's eigenvalues
	particleMPC *pMPC;			//Temporary pointer to MPC particles
	double POPinv,fDIM,invDconst,invBinVol;
	int binPOP,cL;

	avS=0.;
	avS2=0.;
	avS4=0.;
	fDIM = (double) DIM;
	invDconst=1./(fDIM-1.);
	S = malloc ( DIM * sizeof( *S ) );
	for( d=0; d<DIM; d++ ) S[d] = malloc ( DIM * sizeof( *S[d] ) );
	for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) S[i][j] = 0.0;

	for( d=0; d<DIM; d++ ) nBins[d] = XYZ[d]/L;
	if( DIM==_2D) nBins[2]=1;

	for( a=0; a<nBins[0]; a++ ) for( b=0; b<nBins[1]; b++ ) for( c=0; c<nBins[2]; c++ ) {
		//Zero the order parameter tensor
		for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) S[i][j] = 0.;
		binPOP=0;
		//Calculate the order parameter for this Binder bin
		if( DIM==_3D ) cL=(c+1)*L;
		else cL=1;
		for( i=a*L; i<(a+1)*L; i++ ) for( j=b*L; j<(b+1)*L; j++ )  for( k=c*L; k<cL; k++ ) {
			// Find the tensor order parameter for this MPCD cell
			//Calculate the order parameter tensor
			if( CL[i][j][k].pp!=NULL ) {
				binPOP+=CL[i][j][k].POPSRD;
				pMPC = CL[i][j][k].pp;
				//Calculate the order parameter tensor in this cell
				if( LC!=ISOF ) addToTensOrderParam( pMPC,S );
				else addToTensOrderParamVel( pMPC,S );
			}
		}
		POPinv=1./((double) binPOP);
		// Make sum into average
		for( i=0; i<DIM; i++ ) for( j=i; j<DIM; j++ ) S[i][j] *= invDconst*fDIM*POPinv;
		for( i=0; i<DIM; i++ ) S[i][i] -= invDconst;
		// Copy other half
		for( i=0; i<DIM; i++ ) for( j=i+1; j<DIM; j++ ) S[j][i] = S[i][j];
		// From the tensor order parameter find eigenvalues and vectors --- S is written over as normalized eigenvectors
		solveEigensystem( S,DIM,eigval );
		//The scalar order parameter is the largest eigenvalue
		//The scalar order parameter is the largest eigenvalue which is given first ie eigval[0]
		// But can be better approximated (cuz fluctuates about zero) by -2* either of the negative ones (or the average)
		if(DIM==_3D) thisS=-1.*(eigval[1]+eigval[2]);
		else thisS=eigval[0];
		//Add them to averages/sums
		avS+=thisS;
		avS2+=(thisS*thisS);
		avS4+=(thisS*thisS*thisS*thisS);
	}
	invBinVol=1./((double)nBins[0]*nBins[1]*nBins[2]);
	avS*=invBinVol;
	avS2*=invBinVol;
	avS4*=invBinVol;

	UL=1.-avS4/(3.*avS2*avS2);
	return UL;
}

///
/// @brief Apply the BC transformation of the orientation to a MPCD particle caused by a boundary (wall).
///
/// The boundary is considered a ghost particle that affects the orientation of the MPCD particles.
/// If the boundary is movable, there is also a torque between the MPCD and the boundary that is taken care off by applying the appropriate force on the BC due to
/// the reorientation of the MPCD particle.
/// This force would be the opposite torque that the boundary would have applied to the MPCD necesarry for MPCD reorientation.
///
/// @param pp Pointer to an MPCD particle.
/// @param SP List of species.
/// @param WALL Pointer to the considered boundary.
/// @param n The normal of the boundary.
///
void oriBC( particleMPC *pp,spec *SP,bc *WALL,double n[] ) {
	double UN[_3D],UT[_3D],U0[_3D];
	double angleAnch;
	double torque[_3D],r[_3D];			//Torque on MPCD particle --- Not REALLY torque ( angular impulse but time step falls out)
	int i;

	//Zero
	for( i=0; i<_3D; i++ ) {
		U0[i]=0.0;
		UN[i]=0.0;
		UT[i]=0.0;
		r[i]=0.0;
	}

	// Make sure U is a unit vector
	norm(pp->U,DIM);

	//Make U point out of BC
	if(dotprod(pp->U,n,DIM)>=0.0) {
		for(i=0; i<DIM; i++){pp->U[i]*=1.;}
	}
	else{
		for(i=0; i<DIM; i++){pp->U[i]*=-1.;}
	}

	//Calculate the new orientation
	proj( pp->U,n,UN,DIM );		//Calculate normal component of the orientation
	tang( pp->U,UN,UT,DIM );	//Calculate tangential component of orientation
	//Save the original orientation if the BC can move
	if( WALL->DSPLC ) for( i=0; i<DIM; i++ ) {
		U0[i] = pp->U[i];
	}

	//Transform the orientation
	if( feq(WALL->MUN,0.0) && feq(WALL->MUT,0.0) ) {
		//Addition in cartesian coordinates
		for( i=0; i<DIM; i++ ) pp->U[i] = WALL->DUxyz[i];
	}
	else {
		//Multiplication wrt surface normal
		for( i=0; i<DIM; i++ ) {
			UN[i] *= WALL->MUN;
			UT[i] *= WALL->MUT;
		}
		//Combine normal and tangential components
		for( i=0; i<DIM; i++ ) pp->U[i] = UN[i] + UT[i];
		//For measuring K_bend in a pure bend geometry, we want to surpress the z-hat orientation
		pp->U[0] *= WALL->MUxyz[0];
		if( DIM>=_2D ) pp->U[1] *= WALL->MUxyz[1];
		if( DIM>=_3D ) pp->U[2] *= WALL->MUxyz[2];
	}

	// Make sure U0 and U are unit vectors
	norm(U0, _3D);
	norm(pp->U, DIM);

	#ifdef DBG
		if( DBUG == DBGBCORI ) {
			printf( "\tNew u=" );
			pvec(pp->U,DIM);
		}
	#endif

	//Apply torque to BC
	if( WALL->DSPLC ) {
		double u0dotu;
		//Zero
		for( i=0; i<_3D; i++ ) {
			torque[i]=0.0;
		}

		// Calculate angle between initial and final orientation (angleAnch)
		u0dotu = dotprod(U0, pp->U, _3D);
		if (u0dotu >=1.0){
			angleAnch = 0.0; //colloids are crashing without this
		}
		if (u0dotu < 1.0){
			angleAnch = acos(u0dotu); // making sure u0dotu is less than one or the angle blows up
		}

		// Calculating torque on MPCD particle
		crossprod(U0, pp->U, torque);
		norm(torque, _3D); // getting just the direction of torque here

		for (i=0; i<_3D; i++){
			// Multiply unit vector by magnitude.  Now this is torque = (rotational friction coefficient) * (angular velocity). Note: anglar velocity is just the angle.
			torque[i] *= angleAnch;
			torque[i] *= ((SP+pp->SPID)->RFC);
		}

		// Vector between collision point and CM
		for( i=0; i<DIM; i++ ) r[i] = -WALL->Q[i] + pp->Q[i];


		#ifdef DBG
			if( DBUG == DBGBCORI ) {
				printf( "\tTorque on MPC=" );
				pvec(torque,_3D);
			}
		#endif
		//Apply this torque to the BC
		torqueLCBC( WALL, n, U0, torque, (SP+pp->SPID)->LEN, r);
	}
}

///
/// @brief Apply a force on a moveable BC due to orientation boundaries on the MPCD.
///
/// This subroutine applies the appropriate torque on a moveable BC .
/// This subroutine uses the implied torque on the MPCD particle which was necesarry for its reorientation due to the BC boundary conditions.
/// The equal and opposite force to that is the force exerted by MPCD particle on the BC.
/// The force exerted by the BC is converted into rotational and translational motion of the boundary.
///
/// @param WALL pointer to the considered boundary.
/// @param n The normal of the boundary.
/// @param U0 The original orientation of the MPCD particle.
/// @param torqueMPC The torque on the BC caused by the BC exerting a force to reorientate the MPCD.
/// @param rodlength Effective rod length to couple torque on MPC into force on BC (smaller=>stronger; bigger=>weaker) --- json 'len'.
/// @param posColl The vector from the centre of the colloid to the force on the boundary.
///
void torqueLCBC( bc *WALL,double n[], double U0[], double torqueMPC[],double rodlength,  double posColl[] ) {
	int i,j;
	double t_hat[_3D],f_hat[_3D]; //unit vectors: tangent vector to the surface of the colloid (for the tangential component of the force created by the colloid), unit vector of the MPC force acting into the colloid.
	double MPC_force[_3D], force_n[_3D], force_t[_3D]; // MPC force, normal component of the colloid force, tangent component of the colloid force, colloid force.
	double IIwall[_3D][_3D],dL[_3D];									//Inverse of the moment of inertia tensor and angular momentum on BC
	double magT, forceN, forceT;							// magnitude of the MPC torque and magntidude of the normal and tangential components of the colloid force.
	double torqueCol[_3D];										//Torque on colloid (about the centre of rotation/mass of the colloid).

	for( i=0; i<_3D; i++ ) {
		MPC_force[i]=0.0;
		t_hat[i]=0.0;
		f_hat[i]=0.0;
		force_n[i]=0.0;
		force_t[i]=0.0;
		torqueCol[i]=0.0;
		dL[i]=0.0;
		U0[i]=-1.*U0[i]; //flipped it (for preference) as the torques are set up to consider the force pointing towards the boundary.
	}

	norm( U0, _3D); // sometimes this deviates from being a unit vector

	// Define unit vectors
	// 1. Tangent vector (on the colloid surface): If the MPCD particle rotates clockwise (anticlockwise), this vector points clockwise (anticlockwise) around the colloid.

	crossprod( torqueMPC,n,t_hat );
	norm( t_hat, DIM ); //tangent vector at contact point on colloid

	// 3. Direction of force that the colloid provides to the MPCD (to make it rotate)

	crossprod( torqueMPC, U0, f_hat);
	norm( f_hat,_3D);

	// 2. MPCD particle anchoring torque: split into magnitude and unit vector.
	magT = length( torqueMPC,_3D ); // magnitude
	norm( torqueMPC, _3D); // now a unit vector


	#ifdef DBG
		if( DBUG == DBGBCORI ) {
			printf( "\tfhat: " );
			pvec( f_hat,_3D );
		}
	#endif

	// Find magnitude of MPCD rod force: f= (tau/r)*(tau_hat cross r_hat) where r_hat is U0 r is rodlength/2 and tau_hat is the unit vector torqueMPC.
	for(i=0; i<_3D; i++){
		MPC_force[i] = 2.*magT*f_hat[i]/rodlength;
	}

	// Find normal and tangential components of the force given to the colloid by the virtual rod (mirror of MPC particle about surface normal)
	// Force magnitudes
	forceN=dotprod(MPC_force, n, _3D);
	forceT=dotprod(MPC_force, t_hat, _3D);

	// Force vectors
	for( i=0; i<_3D; i++ ){
		force_n[i] = forceN*n[i];
		force_t[i] = -1.*forceT*t_hat[i]; // see also forceT. Without -1, we have found the tangential component of the force on MPC particle, but this tangential component is in the opposite direction to the colloid tangential force.
	}

	#ifdef DBG
		if( DBUG == DBGBCORI ) {
			printf( "\tNormal force on the BC due to rotating an MPCD nematogen: " );
			pvec( force_n,_3D );
			printf( "\tTangential force on the BC due to rotating an MPCD nematogen" );
			pvec( force_t,_3D);
		}
	#endif

	// Translation of colloid: due to normal component of the colloid force
	for( i=0; i<_3D; i++ ) WALL->dV[i] += force_n[i]/(WALL->MASS);

	// Rotation of colloid: about the centre of the colloid due to the tangential component of the colloid force.
	// 1. Calculate the colloid torque about the centre of the colloid
	crossprod( posColl, force_t, torqueCol); //posColl is the vector from the centre of the colloid to the force on the boundary
	// 2. Find the inverse of BC's moment of inertia tensor
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) IIwall[i][j] = 0.0;
	invert3x3(IIwall,WALL->I);
	// 3. Angular impulse from torque and moment of inertia
	dotprodMatVec( IIwall,torqueCol,dL,_3D );
	for( i=0; i<_3D; i++) WALL->dL[i] += dL[i];

	#ifdef DBG
		if( DBUG == DBGBCORI ) {
			printf( "\tResulting Colloid torque:" );
			pvec( torqueCol,_3D );
		}
	#endif
}

/// @brief Perform an MPC collision which conserved angular momentum.
///
/// MPC collision that conserves angular momentum (uses andersen thermostat), and returns the CoM velocity and the local temperature of the cell.
/// See https://pubs.rsc.org/en/content/articlelanding/2015/sm/c5sm00839e
/// Notice that all of this must be entirely in 3D even if system is 2D since angular momentum is a cross product.
/// It takes into account the total mass, velocity and angular momentum of MPC particles, MD particles and swimmers.
///
/// @param CL Class containing cell data (pointed to local cell considered).
/// @param SP List of species.
/// @param SS Pointer to the swimmers.
/// @param KBT Temperature.
/// @param dt Integration timestep.
/// @param CLQ Position of the cell.
/// @param outP Parameter that determines if we should output pressure (0 not, 1 yes).
///
void andersenROT_LC( cell *CL,spec *SP,specSwimmer SS,double KBT,double dt,double *CLQ,int outP ) {
	int i,j,id;
	double MASS;			//MASS
	double RV[CL->POP][_3D];	//Random velocities
	double RS[_3D];			//Sum of random velocities
	double relQ[CL->POP][_3D];	//Relative position
	double diffV[_3D];		//Difference in velocity
	double L[_3D],Llm[_3D],Lnm[_3D];//Angular momentum
	double angterm[_3D];
	double angmom[_3D];
	double W[_3D];
	double VCM[_3D];
	double II[_3D][_3D];		//Inverse of moment of inertia tensor (3D)
	double dp[CL->POP][DIM],relQP[CL->POP][DIM];		//For pressure
	particleMPC *tmpc;		//Temporary particleMPC
	particleMD *tmd;		//Temporary particleMD
	smono *tsm;				//Temporary swimmer monomer

	for( i=0;i<CL->POP;i++ ) for( j=0;j<_3D;j++ ) {
		RV[i][j] = 0.;
		relQ[i][j] = 0.;
	}
	for( j=0;j<_3D;j++ ) {
		RS[j]=0.;
		diffV[j]=0.;
		L[j]=0.;
		angterm[j]=0.;
		angmom[j]=0.;
		W[j]=0.;
		VCM[j] = CL->VCM[j];
		for( i=0;i<_3D;i++ ) II[j][i]=0.;
	}

	/* ****************************************** */
	/* ******* Generate random velocities ******* */
	/* ****************************************** */
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQP[i],dp[i],tmpc,CLQ );
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,MASS );
		for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		MASS = tmd->mass;
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,MASS );
		for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
		tmd = tmd->nextSRD;
		i++;
	}
	//Swimmer particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) MASS = SS.middM;
		else MASS = SS.headM;
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,MASS );
		for( j=0; j<DIM; j++ ) RS[j] += MASS*RV[i][j];
		tsm = tsm->next;
		i++;
	}
	for( j=0; j<DIM; j++ ) RS[j] /= (CL->MASS);

	/* ****************************************** */
	/* ******** Invert moment of inertia ******** */
	/* ****************************************** */
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) II[i][j] = 0.;
	if( DIM == _3D ) {
		// When the population is 2 (don't even get this far if pop=1) then the moment of inertia tensor is singular
		// In that case, can't invert. Therefore, we must neglect angular momentum conservation in such cell
		// This just makes angterm = [0,0,0]
		if( CL->POP < 3 ) for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) II[i][j] = 0.;
		else invert3x3(II,CL->I);
	}
	else if( DIM == _2D ) {
		II[2][2] = 1./CL->I[2][2];
		// The rest don't matter and can remain zero
	}
	else if( DIM == 1 ) {
		//Angular momentum conservation in 1D is nonsequitur. Leave as zero
		II[2][2] = 0.0;
	}
	else {
		printf( "Error: Dimension higher than 3. Change dimension.\n" );
		exit( 1 );
	}

	/* ****************************************** */
	/* ******* Find angular velocity term ******* */
	/* ****************************************** */
	if(DIM>1) {
		for( i=0; i<_3D; i++ ) {
			Llm[i] = 0.;
			Lnm[i] = 0.;
		}
		i=0;
		//MPC particles
		tmpc = CL->pp;
		while( tmpc!=NULL ) {
			//First account for the change in angular momentum due to the effect of the linear momenum collision
			id = tmpc->SPID;
			MASS = (SP+id)->MASS;
			//Position relative to centre of mass
			for( j=0; j<DIM; j++ ) relQ[i][j] = tmpc->Q[j] - CL->CM[j];
			for( j=0; j<DIM; j++ ) diffV[j] = MASS * (tmpc->V[j] - RV[i][j]);
			if( DIM < _3D ) {
				relQ[i][2] = 0.;
				diffV[2] = 0.;
			}
			if( DIM < _2D ) {
				relQ[i][1] = 0.;
				diffV[1] = 0.;
			}
			crossprod( relQ[i],diffV,angmom );
			for( j=0; j<_3D; j++ ) Llm[j] += angmom[j];
			//Next account for the torque on the fluid due to the rotation of the rods
			// net torque=0 (overdamped) therefore HI-torque (on rod) = negative of all other torques
			// torque on fluid = minus HI-torque on rods ---> positive sum of all other torques
			for( j=0; j<_3D; j++ ) Lnm[j] += tmpc->T[j]*dt;

			tmpc = tmpc->next;
			i++;
		}
		//MD particles
		tmd = CL->MDpp;
		while( tmd!=NULL ) {
			//Account for the change in angular momentum due to the effect of the linear momenum collision
			MASS = tmd->mass;
			//Position relative to centre of mass
			relQ[i][0] = tmd->rx - CL->CM[0];
			diffV[0] = MASS * (tmd->vx - RV[i][0]);
			if(DIM > _1D){
				relQ[i][1] = tmd->ry - CL->CM[1];
				diffV[1] = MASS * (tmd->vy - RV[i][1]);
			}
			if( DIM > _2D ){
				relQ[i][2] = tmd->rz - CL->CM[2];
				diffV[2] = MASS * (tmd->vz - RV[i][2]);
			}
			crossprod( relQ[i],diffV,angmom );
			for( j=0; j<_3D; j++ ) Llm[j] += angmom[j];
			tmd = tmd->nextSRD;
			i++;
		}
		//Swimmer particles
		tsm = CL->sp;
		while( tsm!=NULL ) {
			//First account for the change in angular momentum due to the effect of the linear momenum collision
			if( tsm->HorM ) MASS = SS.middM;
			else MASS = SS.headM;
			//Position relative to centre of mass
			for( j=0; j<DIM; j++ ) relQ[i][j] = tsm->Q[j] - CL->CM[j];
			for( j=0; j<DIM; j++ ) diffV[j] = MASS * (tsm->V[j] - RV[i][j]);
			if( DIM < _3D ) {
				relQ[i][2] = 0.;
				diffV[2] = 0.;
			}
			if( DIM < _2D ) {
				relQ[i][1] = 0.;
				diffV[1] = 0.;
			}
			crossprod( relQ[i],diffV,angmom );
			for( j=0; j<_3D; j++ ) Llm[j] += angmom[j];
			tsm = tsm->next;
			i++;
		}
		//Lnm is due to the torque on the rotating rod --- the same must be added to the fluid
		for( i=0; i<_3D; i++ ) L[i]=Llm[i] + Lnm[i];
		//for( i=0; i<_3D; i++ ) L[i]=Llm[i] - Lnm[i];
		dotprodMatVec( II,L,W,_3D );
	}
	else for( j=0; j<DIM; j++ ) W[j]=0.0; // fallback for 1D

	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		MASS = (SP+id)->MASS;
		crossprod( W,relQ[i],angterm );
		for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + RV[i][j] - RS[j] + angterm[j];
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQP[i],dp[i],MASS,tmpc->V,CL );
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		crossprod( W,relQ[i],angterm );
		tmd->vx = VCM[0] + RV[i][0] - RS[0] + angterm[0];
		if( DIM > _1D ){
			tmd->vy = VCM[1] + RV[i][1] - RS[1] + angterm[1];
		}
		if( DIM > _2D ){
			tmd->vz = VCM[2] + RV[i][2] - RS[2] + angterm[2];
		}
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
	//Swimmer particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		crossprod( W,relQ[i],angterm );
		for( j=0; j<DIM; j++ ) tsm->V[j] = VCM[j] + RV[i][j] - RS[j] + angterm[j];
		//Increment link in list
		tsm = tsm->next;
		i++;
	}
}

/// @brief Perform an MPC collision which conserved angular momentum and adds an effective dipole force to the local cell.
///
/// MPC collision that conserves angular momentum (uses andersen thermostat), and returns the CoM velocity and the local temperature of the cell.
/// See https://pubs.rsc.org/en/content/articlelanding/2015/sm/c5sm00839e.
/// Notice that all of this must be entirely in 3D even if system is 2D since angular momentum is a cross product.
/// It takes into account the total mass, velocity and angular momentum of MPC particles, MD particles and swimmers.
/// It also splits the cell in two along the intersection of the centre of mass and perpendicular to the local director.
/// It then adds a velocity component to the particles in opposite direction on either side of this plane to mimic a dipole.
///
/// @param CL Class containing cell data (pointed to local cell considered).
/// @param SP Pointer to the species helper "class" (contains information of the mass, etc, of the considered species).
/// @param SS Pointer to the swimmers.
/// @param KBT Temperature.
/// @param RELAX At the moment doesn't do anything.
/// @param RTECH Input which defines different rotation techniques.
/// @param dt Integration timestep.
/// @param CLQ Position of the cell.
/// @param outP Parameter that determines if we should output pressure (0 not, 1 yes).
///
void dipoleAndersenROT_LC( cell *CL,spec *SP,specSwimmer SS,double KBT,double RELAX,double dt,int RTECH,double *CLQ,int outP ) {
	int i,j,id;
	double M,MASS,ACT,pmOne,sigWidth,sigPos;
	double RV[CL->POP][_3D];	//Random velocities
	double RS[_3D];			//Sum of random velocities
	double DV[CL->POP][_3D];	//Damping velocities
	double AV[CL->POP][DIM];	//Active velocities
	double AS[DIM];			//Sum of active velocities
	double relQ[CL->POP][_3D];	//Relative position
	double diffV[_3D];		//Difference in velocity
	double L[_3D],Llm[_3D],Lnm[_3D];//Angular momentum
	double angterm[_3D];
	double angmom[_3D];
	double W[_3D];
	double VCM[_3D];
	double II[_3D][_3D];		//Inverse of moment of inertia tensor (3D)
	double dp[CL->POP][DIM],relQP[CL->POP][DIM];		//For pressure
	double pW;			//The particle's pW for passing the plane
	bc PLANE;			//The plane that cuts the cell in half
	particleMPC *tmpc;		//Temporary particleMPC
	particleMD *tmd;		//Temporary particleMD
	smono *tsm;			//Temporary swimmer mnomer

	// Zero arrays
	for( i=0;i<CL->POP;i++ ) for( j=0;j<_3D;j++ ) {
		RV[i][j] = 0.;
		DV[i][j] = 0.;
		AV[i][j] = 0.;
		relQ[i][j] = 0.;
	}
	for( j=0;j<_3D;j++ ) {
		RS[j]=0.;
		AS[j]=0.;
		diffV[j]=0.;
		L[j]=0.;
		angterm[j]=0.;
		angmom[j]=0.;
		W[j]=0.;
		VCM[j] = CL->VCM[j];
	}
	for( i=0;i<_3D;i++ ) for( j=0;j<_3D;j++ ) II[i][j]=0.;

	//Define the plane normal to the centre of mass velocity at the centre of mass position
	for( i=0;i<4;i++ ) PLANE.P[i]=1;
	PLANE.INV=0;
	PLANE.ABS=0;
	PLANE.R=0.0;
	PLANE.ROTSYMM[0]=4.0;
	PLANE.ROTSYMM[1]=4.0;
	//Normal
	for( i=0; i<DIM; i++ ) PLANE.A[i] = CL->DIR[i];
	//Position
	for( i=0; i<DIM; i++ ) PLANE.Q[i] = CL->CM[i];

	MASS=CL->MASS;
	//Calculate total activity and the average sigmoidal width of cell
	tmpc = CL->pp;
	ACT=0.;
	sigWidth=0.;
	sigPos=0.;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;

		// do a check to see if the subpopulation is of sufficient quantity to compute activity
		if ( ((double)(SP+id)->MINACTRATIO == 0.0) || (((double)(SP+id)->MINACTRATIO)*nDNST < (double)CL->SP[id]) ) {
			ACT += (double)(SP+id)->ACT / (double)(SP+id)->MASS;
			sigWidth += (double)(SP+id)->SIGWIDTH;
			sigPos += (double)(SP+id)->SIGPOS;
		}
		tmpc = tmpc->next;
	}
	//If DIPOLE_DIR_SUM then use the sum just calculated
	//If DIPOLE_DIR_AV or DIPOLE_DIR_SIG then use the average value everywhere
	if( RTECH==DIPOLE_DIR_AV || RTECH==DIPOLE_DIR_SIG) ACT *= nDNST/((double)CL->POP);

	// Now, if using sigmoidal dipole set up a sigmoidal falloff based on the cell population
	if (RTECH==DIPOLE_DIR_SIG || RTECH==DIPOLE_DIR_SIG_SUM) {
		// compute average within cell by normalising with cell population
		sigWidth /= (double)CL->POP;
		sigPos /= (double)CL->POP;

		// compute the sigmoidal falloff function
		double falloffFactor = (1 - tanh( ((double)CL->POP  - nDNST * (1 + sigPos) ) / (nDNST * sigWidth) ) );
		double rescaleFactor = 2;

		// rescale the activity
		ACT *= falloffFactor / rescaleFactor;
	}

	// Scale activity by the timestep size to remove timestep dependence
	ACT *= dt;

	/* ****************************************** */
	/* ******* Generate random velocities ******* */
	/* ****************************************** */
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		M = (SP+id)->MASS;
		//Pressure term
		if( outP ) calcPressureColl_preColl( relQP[i],dp[i],tmpc,CLQ );
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,M );
		for( j=0; j<DIM; j++ ) RS[j] += M*RV[i][j];
		for( j=0; j<DIM; j++ ) DV[i][j] = ((SP+id)->DAMP)*(CL->VCM[j]);
		if( fneq(ACT,0.0) ) {
			//Check which side of the plane
			pW = calcW( PLANE,*tmpc );
			if( pW<=0 ) pmOne=-1.;
			else pmOne=1.;
			for( j=0; j<DIM; j++ ) AV[i][j] = pmOne*ACT*PLANE.A[j];
			for( j=0; j<DIM; j++ ) AS[j] += M*AV[i][j];
		}
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		M = tmd->mass;
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,M );
		for( j=0; j<DIM; j++ ) RS[j] += M*RV[i][j];
		tmd = tmd->nextSRD;
		i++;
	}
	//Swimmer monomer
	tsm = CL->sp;
	while( tsm!=NULL ) {
		if( tsm->HorM ) M = SS.middM;
		else M = SS.headM;
		for( j=0; j<DIM; j++ ) RV[i][j] = genrand_gaussMB( KBT,M );
		for( j=0; j<DIM; j++ ) RS[j] += M*RV[i][j];
		if( fneq(ACT,0.0) ) {
			//Check which side of the plane
			pW = calcW( PLANE,*tmpc );
			if( pW<=0 ) pmOne=-1.;
			else pmOne=1.;
			for( j=0; j<DIM; j++ ) AV[i][j] = pmOne*ACT*PLANE.A[j];
			for( j=0; j<DIM; j++ ) AS[j] += M*AV[i][j];
		}
		tsm = tsm->next;
		i++;
	}
	//Normalize
	for( j=0; j<DIM; j++ ) RS[j] /= MASS;
	for( j=0; j<DIM; j++ ) AS[j] /= MASS;

	/* ****************************************** */
	/* ******** Invert moment of inertia ******** */
	/* ****************************************** */
	for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) II[i][j] = 0.;
	if( DIM == _3D ) {
		// When the population is 2 (don't even get this far if pop=1) then the moment of inertia tensor is singular
		// In that case, can't invert. Therefore, we must neglect angular momentum conservation in such cell
		// This just makes angterm = [0,0,0]
		if( CL->POP < 3 ) for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) II[i][j] = 0.;
		else invert3x3(II,CL->I);
	}
	else if( DIM == _2D ) {
		II[2][2] = 1./CL->I[2][2];
		// The rest don't matter and can remain zero
	}
	else if( DIM == 1 ) {
		//Angular momentum conservation in 1D is nonsequitur. Leave as zero
		II[2][2] = 0.0;
	}
	else {
		printf( "Error: Dimension higher than 3. Change dimension.\n" );
		exit( 1 );
	}

	/* ****************************************** */
	/* ******* Find angular velocity term ******* */
	/* ****************************************** */
	if(DIM>1) {
		for( i=0; i<_3D; i++ ) {
			Llm[i] = 0.;
			Lnm[i] = 0.;
		}
		i=0;
		//MPC particles
		tmpc = CL->pp;
		while( tmpc!=NULL ) {
			//First account for the change in angular momentum due to the effect of the linear momenum collision
			id = tmpc->SPID;
			M = (SP+id)->MASS;
			//Position relative to centre of mass
			for( j=0; j<DIM; j++ ) relQ[i][j] = tmpc->Q[j] - CL->CM[j];
			//Difference in velocity --- with BOTH active and random terms
			for( j=0; j<DIM; j++ ) diffV[j] = M * (tmpc->V[j] - RV[i][j] - AV[i][j]);
			if( DIM < _3D ) {
				relQ[i][2] = 0.;
				diffV[2] = 0.;
			}
			if( DIM < _2D ) {
				relQ[i][1] = 0.;
				diffV[1] = 0.;
			}
			crossprod( relQ[i],diffV,angmom );
			for( j=0; j<_3D; j++ ) Llm[j] += angmom[j];
			//Next account for the torque on the fluid due to the rotation of the rods
			// net torque=0 (overdamped) therefore HI-torque (on rod) = negative of all other torques
			// torque on fluid = minus HI-torque on rods ---> positive sum of all other torques
			for( j=0; j<_3D; j++ ) Lnm[j] += tmpc->T[j]*dt;

			tmpc = tmpc->next;
			i++;
		}
		//MD particles
		tmd = CL->MDpp;
		while( tmd!=NULL ) {
			//Account for the change in angular momentum due to the effect of the linear momenum collision
			M = tmd->mass;
			//Position relative to centre of mass
			relQ[i][0] = tmd->rx - CL->CM[0];
			diffV[0] = MASS * (tmd->vx - RV[i][0]);
			if(DIM > _1D){
				relQ[i][1] = tmd->ry - CL->CM[1];
				diffV[1] = MASS * (tmd->vy - RV[i][1]);
			}
			if( DIM > _2D ){
				relQ[i][2] = tmd->rz - CL->CM[2];
				diffV[2] = MASS * (tmd->vz - RV[i][2]);
			}
			crossprod( relQ[i],diffV,angmom );
			for( j=0; j<_3D; j++ ) Llm[j] += angmom[j];
			tmd = tmd->nextSRD;
			i++;
		}
		//Swimmer particles
		tsm = CL->sp;
		while( tsm!=NULL ) {
			//First account for the change in angular momentum due to the effect of the linear momenum collision
			if( tsm->HorM ) M = SS.middM;
			else M = SS.headM;
			//Position relative to centre of mass
			for( j=0; j<DIM; j++ ) relQ[i][j] = tsm->Q[j] - CL->CM[j];
			//Difference in velocity
			for( j=0; j<DIM; j++ ) diffV[j] = M * (tsm->V[j] - RV[i][j] );
			if( DIM < _3D ) {
				relQ[i][2] = 0.;
				diffV[2] = 0.;
			}
			if( DIM < _2D ) {
				relQ[i][1] = 0.;
				diffV[1] = 0.;
			}
			crossprod( relQ[i],diffV,angmom );
			for( j=0; j<_3D; j++ ) Llm[j] += angmom[j];
			tsm = tsm->next;
			i++;
		}

		for( i=0; i<_3D; i++ ) L[i]=Llm[i]+Lnm[i];
		dotprodMatVec( II,L,W,_3D );
	}
	else for( j=0; j<DIM; j++ ) W[j]=0.0; // fallback for 1D

	/* ****************************************** */
	/* *************** Collision **************** */
	/* ****************************************** */
	i=0;
	//MPC particles
	tmpc = CL->pp;
	while( tmpc!=NULL ) {
		id = tmpc->SPID;
		M = (SP+id)->MASS;
		crossprod( W,relQ[i],angterm );
		for( j=0; j<DIM; j++ ) tmpc->V[j] = VCM[j] + RV[i][j] - RS[j] + angterm[j]  + AV[i][j] - AS[j] -DV[i][j];
		//Pressure term
		if( outP ) calcPressureColl_postColl( relQP[i],dp[i],M,tmpc->V,CL );
		//Increment link in list
		tmpc = tmpc->next;
		i++;
	}
	//MD particles
	tmd = CL->MDpp;
	while( tmd!=NULL ) {
		crossprod( W,relQ[i],angterm );
		tmd->vx = VCM[0] + RV[i][0] - RS[0] + angterm[0] ;
		if( DIM > _1D ){
			tmd->vy = VCM[1] + RV[i][1] - RS[1] + angterm[1];
		}
		if( DIM > _2D ){
			tmd->vz = VCM[2] + RV[i][2] - RS[2] + angterm[2];
		}
		//Increment link in list
		tmd = tmd->nextSRD;
		i++;
	}
	//Swimmer particles
	tsm = CL->sp;
	while( tsm!=NULL ) {
		crossprod( W,relQ[i],angterm );
		for( j=0; j<DIM; j++ ) tsm->V[j] = VCM[j] + RV[i][j] - RS[j] + angterm[j]  + AV[i][j] - AS[j] -DV[i][j];
		//Increment link in list
		tsm = tsm->next;
		i++;
	}
}

///
/// @brief Find the local velocity gradient tensor of all cells.
///
/// This routine finds the local velocity gradient tensor E[i][j]= dv[i]/dx[j] of all cells.
/// It does not bother with periodic BCs.
/// It uses different routines for the different dimensions for convenience (currently used the D3Q15 and D2Q9 neighbours. Does not work with 1D).
///
/// @param CL All of the MPCD cells (including the linked list of particles in each cell). 
///
void localVelGrad( cell ***CL ) {
	// if( DIM == _3D ) velGrad3D( CL );
	if( DIM == _3D ) velGradD3Q15( CL );
// 	else if( DIM == _2D ) velGrad2D( CL );
	else if( DIM == _2D ) velGradD2Q9( CL );
// 	else if( DIM == _1D ) velGrad1D( CL );
}

///
/// @brief Find the local velocity gradient tensor of all cells for 3D using D3Q15 neighbour grid.
///
/// This routine finds the local velocity gradient tensor E[i][j]= dv[i]/dx[j] of all cells for 3D using D3Q15 neighbour grid.
/// It does not bother with periodic BCs and doesn't calculate corners.
///
/// @param CL All of the MPCD cells (including the linked list of particles in each cell). 
///
void velGradD3Q15( cell ***CL ) {

	int a,b,c,i,j,k;
	int numNodes = 14;
	int unitVec[14][3] = { {1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1} };		 	//LB nodal unit vectors
	double w123456=1./9., w789=1./72.;
	double w[14]={ w123456,w123456,w123456,w123456,w123456,w123456,w789,w789,w789,w789,w789,w789,w789,w789 };	//LB weighting
	double Tinv=3.;

	//Bulk
	for( a=1; a<XYZ[0]; a++ ) for( b=1; b<XYZ[1]; b++ ) for( c=1; c<XYZ[2]; c++ ) for( i=0; i<DIM; i++ ) {
		//D3Q15 Gradient
		for( j=0; j<DIM; j++ ) {
			CL[a][b][c].E[i][j] = 0.;
			for( k=0; k<numNodes; k++ ) {
				CL[a][b][c].E[i][j] += w[k] * unitVec[k][j] * CL[ a+unitVec[k][0] ][ b+unitVec[k][1] ][ c+unitVec[k][2] ].VCM[i];
			}
			CL[a][b][c].E[i][j] *= Tinv;
		}
	}
	//x-edge
	a=0;
	for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) {
		CL[a][b][c].E[i][j] =CL[a+1][b][c].E[i][j];
	}
	a=XYZ[0];
	for( b=0; b<XYZ[1]; b++ ) for( c=0; c<XYZ[2]; c++ ) for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) {
		CL[a][b][c].E[i][j] =CL[a-1][b][c].E[i][j];
	}
	//y-edge
	b=0;
	for( a=0; a<XYZ[0]; a++ ) for( c=0; c<XYZ[2]; c++ ) for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) {
		CL[a][b][c].E[i][j] =CL[a][b+1][c].E[i][j];
	}
	b=XYZ[1];
	for( a=0; a<XYZ[0]; a++ ) for( c=0; c<XYZ[2]; c++ ) for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) {
		CL[a][b][c].E[i][j] =CL[a][b-1][c].E[i][j];
	}
	//z-edge
	c=0;
	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) {
		CL[a][b][c].E[i][j] =CL[a][b][c+1].E[i][j];
	}
	c=XYZ[2];
	for( a=0; a<XYZ[0]; a++ ) for( b=0; b<XYZ[1]; b++ ) for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) {
		CL[a][b][c].E[i][j] =CL[a][b][c-1].E[i][j];
	}
	//Corners
	//Could do average of neighbours
}

///
/// @brief Find the local velocity gradient tensor of all cells for 2D using D2Q9 neighbour grid.
///
/// This routine finds the local velocity gradient tensor E[i][j]= dv[i]/dx[j] of all cells for 2D using D3Q9 neighbour grid.
/// It does not bother with periodic BCs and doesn't calculate corners.
/// It uses LB-type derivatives to get isotropy correct to leading order - see https://iopscience.iop.org/article/10.1209/0295-5075/101/50006 .
///
/// @param CL All of the MPCD cells (including the linked list of particles in each cell).
///
void velGradD2Q9( cell ***CL ) {
	int a,b,c,i,j,k;
	int numNodes = 8;
	int unitVec[8][3] = { {1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{1,1,0},{-1,1,0},{-1,-1,0},{1,-1,0} }; 	//LB nodal unit vectors
	double w1234=1./9., w5678=1./36.;
	double w[8]={ w1234,w1234,w1234,w1234,w5678,w5678,w5678,w5678 };				//LB weighting
	double Tinv=3.;

	c=0;		//There is no z
	//Bulk
	for( a=1; a<XYZ[0]; a++ ) for( b=1; b<XYZ[1]; b++ ) for( i=0; i<DIM; i++ ) {
		//D2Q9 Gradient
		for( j=0; j<DIM; j++ ) {
			CL[a][b][c].E[i][j] = 0.;
			for( k=0; k<numNodes; k++ ) {
				CL[a][b][c].E[i][j] += w[k] * unitVec[k][j] * CL[ a+unitVec[k][0] ][ b+unitVec[k][1] ][ c+unitVec[k][2] ].VCM[i];
			}
			CL[a][b][c].E[i][j] *= Tinv;
		}
	}
	//x-edge
	a=0;
	for( b=0; b<XYZ[1]; b++ ) for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) CL[a][b][c].E[i][j] =CL[a+1][b][c].E[i][j];
	a=XYZ[0];
	for( b=0; b<XYZ[1]; b++ ) for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) CL[a][b][c].E[i][j] =CL[a-1][b][c].E[i][j];
	//y-edge
	b=0;
	for( a=0; a<XYZ[0]; a++ ) for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) CL[a][b][c].E[i][j] = CL[a][b+1][c].E[i][j];
	b=XYZ[1];
	for( a=0; a<XYZ[0]; a++ ) for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) CL[a][b][c].E[i][j] = CL[a][b-1][c].E[i][j];
	//Corners
	//Could do average of neighbours
}

///
/// @brief Find the local velocity gradient tensor of all cells for 3D
///
/// Currently not in use.
/// This routine finds the local velocity gradient tensor E[i][j]= dv[i]/dx[j] of all cells for 3D using the Carthesian neighbours.
///
/// @param CL All of the MPCD cells (including the linked list of particles in each cell).
///
void velGrad3D( cell ***CL ) {
	int a,b,c,i;
	//Bulk
	for( a=1; a<XYZ[0]; a++ ) for( b=1; b<XYZ[1]; b++ ) for( c=1; c<XYZ[2]; c++ ) for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = centredDeriv( CL[a-1][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = centredDeriv( CL[a][b-1][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = centredDeriv( CL[a][b][c-1].VCM[i],CL[a][b][c+1].VCM[i],1. );
	}
	//x-edge
	a=0;
	for( b=1; b<XYZ[1]; b++ ) for( c=1; c<XYZ[2]; c++ ) for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = forwardDeriv( CL[a][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = centredDeriv( CL[a][b-1][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = centredDeriv( CL[a][b][c-1].VCM[i],CL[a][b][c+1].VCM[i],1. );
	}
	a=XYZ[0];
	for( b=1; b<XYZ[1]; b++ ) for( c=1; c<XYZ[2]; c++ ) for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = backwardDeriv( CL[a][b][c].VCM[i],CL[a-1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = centredDeriv( CL[a][b-1][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = centredDeriv( CL[a][b][c-1].VCM[i],CL[a][b][c+1].VCM[i],1. );
	}
	//y-edge
	b=0;
	for( a=1; a<XYZ[0]; a++ ) for( c=1; c<XYZ[2]; c++ ) for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = centredDeriv( CL[a-1][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = forwardDeriv( CL[a][b][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = centredDeriv( CL[a][b][c-1].VCM[i],CL[a][b][c+1].VCM[i],1. );
	}
	b=XYZ[1];
	for( a=1; a<XYZ[0]; a++ ) for( c=1; c<XYZ[2]; c++ ) for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = centredDeriv( CL[a-1][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = backwardDeriv( CL[a][b][c].VCM[i],CL[a][b-1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = centredDeriv( CL[a][b][c-1].VCM[i],CL[a][b][c+1].VCM[i],1. );
	}
	//z-edge
	c=0;
	for( a=1; a<XYZ[0]; a++ ) for( b=1; b<XYZ[1]; b++ ) for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = centredDeriv( CL[a-1][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = centredDeriv( CL[a][b][c-1].VCM[i],CL[a][b+1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = forwardDeriv( CL[a][b][c].VCM[i],CL[a][b][c+1].VCM[i],1. );
	}
	c=XYZ[2];
	for( a=1; a<XYZ[0]; a++ ) for( b=1; b<XYZ[1]; b++ ) for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = centredDeriv( CL[a-1][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = centredDeriv( CL[a][b-1][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = backwardDeriv( CL[a][b][c].VCM[i],CL[a][b][c-1].VCM[i],1. );
	}
	//Corners
	a=0;
	b=0;
	c=0;			//a=0,b=0,c=0
	for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = forwardDeriv( CL[a][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = forwardDeriv( CL[a][b][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = forwardDeriv( CL[a][b][c].VCM[i],CL[a][b][c+1].VCM[i],1. );
	}
	a=XYZ[0];		//a=X,b=0,c=0
	for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = backwardDeriv( CL[a][b][c].VCM[i],CL[a-1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = forwardDeriv( CL[a][b][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = forwardDeriv( CL[a][b][c].VCM[i],CL[a][b][c+1].VCM[i],1. );
	}
	b=XYZ[1];		//a=X,b=Y,c=0
	for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = backwardDeriv( CL[a][b][c].VCM[i],CL[a-1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = backwardDeriv( CL[a][b][c].VCM[i],CL[a][b-1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = forwardDeriv( CL[a][b][c].VCM[i],CL[a][b][c+1].VCM[i],1. );
	}
	c=XYZ[2];		//a=X,b=Y,c=Z
	for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = backwardDeriv( CL[a][b][c].VCM[i],CL[a-1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = backwardDeriv( CL[a][b][c].VCM[i],CL[a][b-1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = backwardDeriv( CL[a][b][c].VCM[i],CL[a][b][c-1].VCM[i],1. );
	}
	b=0;			//a=X,b=0,c=Z
	for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = backwardDeriv( CL[a][b][c].VCM[i],CL[a-1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = forwardDeriv( CL[a][b][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = backwardDeriv( CL[a][b][c].VCM[i],CL[a][b][c-1].VCM[i],1. );
	}
	a=0;			//a=0,b=0,c=Z
	for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = forwardDeriv( CL[a][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = forwardDeriv( CL[a][b][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = backwardDeriv( CL[a][b][c].VCM[i],CL[a][b][c-1].VCM[i],1. );
	}
	b=XYZ[1];		//a=0,b=Y,c=Z
	for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = forwardDeriv( CL[a][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = backwardDeriv( CL[a][b][c].VCM[i],CL[a][b-1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = backwardDeriv( CL[a][b][c].VCM[i],CL[a][b][c-1].VCM[i],1. );
	}
	c=0;			//a=0,b=Y,c=0
	for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = forwardDeriv( CL[a][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = backwardDeriv( CL[a][b][c].VCM[i],CL[a][b-1][c].VCM[i],1. );
		CL[a][b][c].E[i][2] = forwardDeriv( CL[a][b][c].VCM[i],CL[a][b][c+1].VCM[i],1. );
	}
}

///
/// @brief Find the local velocity gradient tensor of all cells for 2D.
///
/// Currently not in use.
/// This routine finds the local velocity gradient tensor E[i][j]= dv[i]/dx[j] of all cells for 2D using the Carthesian neighbours.
/// There is no z.
///
/// @param CL All of the MPCD cells (including the linked list of particles in each cell).
///
void velGrad2D( cell ***CL ) {
	int a,b,c,i;
	c=0;		//There is no z
	//Bulk
	for( a=1; a<XYZ[0]; a++ ) for( b=1; b<XYZ[1]; b++ ) for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = centredDeriv( CL[a-1][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = centredDeriv( CL[a][b-1][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
	}
	//x-edge
	a=0;
	for( b=1; b<XYZ[1]; b++ ) for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = forwardDeriv( CL[a][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = centredDeriv( CL[a][b-1][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
	}
	a=XYZ[0];
	for( b=1; b<XYZ[1]; b++ ) for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = backwardDeriv( CL[a][b][c].VCM[i],CL[a-1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = centredDeriv( CL[a][b-1][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
	}
	//y-edge
	b=0;
	for( a=1; a<XYZ[0]; a++ ) for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = centredDeriv( CL[a-1][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = forwardDeriv( CL[a][b][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
	}
	b=XYZ[1];
	for( a=1; a<XYZ[0]; a++ ) for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = centredDeriv( CL[a-1][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = backwardDeriv( CL[a][b][c].VCM[i],CL[a][b-1][c].VCM[i],1. );
	}
	//Corners
	a=0;
	b=0;			//a=0,b=0,c=0
	for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = forwardDeriv( CL[a][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = forwardDeriv( CL[a][b][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
	}
	a=XYZ[0];		//a=X,b=0,c=0
	for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = backwardDeriv( CL[a][b][c].VCM[i],CL[a-1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = forwardDeriv( CL[a][b][c].VCM[i],CL[a][b+1][c].VCM[i],1. );
	}
	b=XYZ[1];		//a=X,b=Y,c=0
	for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = backwardDeriv( CL[a][b][c].VCM[i],CL[a-1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = backwardDeriv( CL[a][b][c].VCM[i],CL[a][b-1][c].VCM[i],1. );
	}
	a=0;			//a=0,b=Y,c=0
	for( i=0; i<DIM; i++ ) {
		CL[a][b][c].E[i][0] = forwardDeriv( CL[a][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
		CL[a][b][c].E[i][1] = backwardDeriv( CL[a][b][c].VCM[i],CL[a][b-1][c].VCM[i],1. );
	}
}

///
/// @brief Find the local velocity gradient tensor of all cells for 1D.
///
/// Currently not in use.
/// This routine finds the local velocity gradient tensor E[i][j]= dv[i]/dx[j] of all cells for 1D using the Carthesian neighbours.
///
/// @param CL All of the MPCD cells (including the linked list of particles in each cell).
///
void velGrad1D( cell ***CL ) {

//   There is no z or y --- this is just silly
	int a,b,c,i;
	c=0;		//There is no z
	b=0;
	//Bulk
	for( a=1; a<XYZ[0]; a++ ) for( i=0; i<DIM; i++ ) CL[a][b][c].E[i][0] = centredDeriv( CL[a-1][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
	//x-edge
	a=0;
	for( i=0; i<DIM; i++ ) 	CL[a][b][c].E[i][0] = forwardDeriv( CL[a][b][c].VCM[i],CL[a+1][b][c].VCM[i],1. );
	a=XYZ[0];
	for( i=0; i<DIM; i++ ) CL[a][b][c].E[i][0] = backwardDeriv( CL[a][b][c].VCM[i],CL[a-1][b][c].VCM[i],1. );
}

///
/// @brief Calculate and return the local topological charge.
///
/// A subroutine that finds and gives the topological charge (in the 2D XY-plane) at the local cell 
/// with position (i,j,k) by rotating the cell directory around the 8 neighbours of the cell.
///
/// @param CL All of the MPCD cells (including the linked list of particles in each cell).
/// @param i The x index of cell considered.
/// @param j The y index of cell considered.
/// @param k The z index of cell considered.
/// @return Returns the topological charge at position (i,j,k).
///
double topoChargeLocal( cell ***CL, int i, int j, int k){
	double phi = .0; //reset

	phi += topoSmallestAngle( CL[i+1][j][k].DIR, CL[i+1][j+1][k].DIR);
	phi += topoSmallestAngle( CL[i+1][j+1][k].DIR, CL[i][j+1][k].DIR);
	phi += topoSmallestAngle( CL[i][j+1][k].DIR, CL[i-1][j+1][k].DIR);
	phi += topoSmallestAngle( CL[i-1][j+1][k].DIR, CL[i-1][j][k].DIR);
	phi += topoSmallestAngle( CL[i-1][j][k].DIR, CL[i-1][j-1][k].DIR);
	phi += topoSmallestAngle( CL[i-1][j-1][k].DIR, CL[i][j-1][k].DIR);
	phi += topoSmallestAngle( CL[i][j-1][k].DIR, CL[i+1][j-1][k].DIR);
	phi += topoSmallestAngle( CL[i+1][j-1][k].DIR, CL[i+1][j][k].DIR);

	return 0.5*phi/pi;
}

///
/// @brief Return the smallest vector between two directors u and v.
///
/// Takes into account that the director is defined as u=-u.
///
/// @param u A vector (director).
/// @param v A vector (director).
/// @return Returns the (smallest) angle between the two directors.
///
double topoSmallestAngle( double u[], double v[]){
	//copy pasting code from python here to calculate smallest angle
	double dot = dotprod(u, v, _2D);
	double det = u[0]*v[1] - u[1]*v[0];
	double phi = atan2(fabs(det), dot);

	if (phi > 0.5*pi){
		v[0] = -v[0];
		v[1] = -v[1];
		dot = dotprod(u, v, _2D);
		det = u[0]*v[1] - u[1]*v[0];
	}
	double sign = 1.0;
	if (det < 0) sign = -1.0;

	return sign*atan2(fabs(det), dot);
}

///
/// @brief Return the angle of a topological defect.
///
/// A defect first has to be detected using topoChargeLocal() as the charge of that function needs to be given in this function.
/// This function calculated the derivatives of the Q tensor of the local cell and its neighbours and then uses the equation  taken from: https://pubs.rsc.org/en/content/articlelanding/2016/sm/c6sm01146b/ .
/// Currently doesn't seem to be able to handle when the cell are taken at the boundary.
///
/// @param CL All of the MPCD cells (including the linked list of particles in each cell).
/// @param x The x index of cell considered.
/// @param y The y index of cell considered.
/// @param z The z index of cell considered.
/// @param charge Charge of the defect.
/// @returns Returns the angle of a defect.
/// @see topoChargeLocal()
///
double topoAngleLocal( cell ***CL, int x, int y, int z, double charge){
	//init summing vars
	double sumTop = 0.0;
	double sumBot = 0.0;
	// prepare sign based on charge
	int sign = 1;
	if (charge < 0) sign = -1;

	//loop through neighbouring cells to compute average
	int i, j = 0;
	for( i = x-1; i < x+2; i++) for( j = y-1; j < y+2; j++) if(!( (i==x) && (j==y))){
		//compute necessary partials of this NEIGHBOURING cell
		//we need partials in x and y of Q_{xx} and Q_{xy}, so compute them using finite central diff
		//first, get the necessary Q tensors
		double QTop[_2D][_2D];
		computeQ(CL[i][j+1][z], QTop);
		double QBot[_2D][_2D];
		computeQ(CL[i][j-1][z], QBot);
		double QLeft[_2D][_2D];
		computeQ(CL[i-1][j][z], QLeft);
		double QRight[_2D][_2D];
		computeQ(CL[i+1][j][z], QRight);

		//now compute derivatives of format (partial axis) Q (Q element)
		double xQxx = centredDeriv(QLeft[0][0], QRight[0][0], 1.0);
		double xQxy = centredDeriv(QLeft[1][0], QRight[1][0], 1.0);
		double yQxx = centredDeriv(QBot[0][0], QTop[0][0], 1.0);
		double yQxy = centredDeriv(QBot[1][0], QTop[1][0], 1.0);
		
		sumTop += sign*xQxy - yQxx; // compute top here
		sumBot += xQxx + sign*yQxy; // compute bot here

		/***TODO: test computeQ() against tensOrderParam() above a _high_ tolerance***/
	}

	double angle = (charge / (1.0 - charge)) * atan2(sumTop, sumBot); //compute angle per the equation
	return angle;
}

///
/// @brief Calculate the Q tensor.
///
/// Calculates the Q tensor from the director and nematic order in CL.
/// Currently only works in 2D.
///
/// @param CL Contains cell data at previously defined indices.
/// @param output The Q tensor.
///
void computeQ(cell CL, double output[_2D][_2D]){
	//FIXME: only works for 2D for now!!!
	// set up the Q tensor object we plan to return
	int i, j = 0;
	for( i = 0; i < DIM; i++) for( j = 0; j < DIM; j++) output[i][j] = 0.0;

	//using personal outer product here because of issues with casting 2D arrays
	for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) output[i][j] = CL.DIR[i]*CL.DIR[j];
	for( i = 0; i < DIM; i++) output[i][i] -= 1.0/((double)DIM); //subtract I term
	for( i = 0; i < DIM; i++) for( j = 0; j < DIM; j++) output[i][j] *= CL.S; // scale with order parameter
}
