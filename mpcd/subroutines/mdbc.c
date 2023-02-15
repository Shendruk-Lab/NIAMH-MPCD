///
/// @file
///
/// @brief Boundary conditions (BCs) on molecular dynamics (MD) particles.
///
/// This file includes the routines that determine the interaction of the MD particles with
/// the BCs. For instance how position and velocity of the MD particles are updated when  
/// crossing a bounday and if the energy, momentum and angular momentum are conserved.
/// 

# include<stdio.h>
# include<math.h>
# include<time.h>
# include<string.h>
# include<stdlib.h>

# include "../headers/definitions.h"
# include "../headers/globals.h"
# include "../headers/SRDclss.h"
# include "../headers/rand.h"
# include "../headers/pout.h"
# include "../headers/ctools.h"
# include "../headers/mtools.h"
# include "../headers/pout.h"
# include "../headers/mpc.h"
# include "../headers/therm.h"
# include "../headers/bc.h"
# include "../headers/init.h"
# include "../headers/mdbc.h"

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ******** MD PARTICLES PASSING BCs ******** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

///
/// @brief	Applies the BCs to the MD particle, in case it is required.
///
/// It checks if the particle is inside the boundaries. If it's inside the boundaries it does not
/// do anything.
/// But if it is not, it rewinds it back to its old position, then it calculates the
/// time taken for the particle to collide with the boundary all through chooseBC_MD().
/// It streams that time, collides with the boundary, then it streams in to the control volume for the 
/// rest of the streaming time.
/// If the function fails to bring it inside, it is simply brought back to its old position with no
///	velocity.
///
/// @param atom		The MD particle.
/// @param WALL		All of the walls (BCs) that the particle might interact with. 
/// @param KBT		Thermal energy.
/// @param t_step	The MD timestep increment.
/// @see 			MPC_BCcollision()
/// 
void MD_BCcollision( particleMD *atom,bc WALL[],double KBT,double t_step ) {
	double t_delta;			//time passed so far
	double time;			//time left to move for
	double t_min=0.;		//smallest time
	int chosenBC=0;			//Particle to go with t_min
	int flag = 1;			//flag for if should keep looping. Loop while 1.
	int cnt = 0;
	double n[_3D] = {0.,0.,0.};	//Normal to the surface
	double W, W1 = 0.0;
	double shift[DIM];
	double RX=0.0,RY=0.0,RZ=0.0;
	double WX=0.0,WY=0.0,WZ=0.0;

	t_delta = 0.;
	time = t_step;

	while( flag ) {
		//We must check if the particle is inside any of the BCs
		chooseBC_MD( WALL,atom,&t_min,&W,&chosenBC,time,t_step );

		//If no particles were inside then we are done.
		if( W > -TOL ) flag = 0;
		//Otherwise, COLLISON
		else {
			cnt++;
			//We have the BC to collide with and the time at which it collided
			//Rewind the particle back to it's old position
			rewind_MD( atom,time );
			if(cnt==1){
				RX=atom->rx;
				RY=atom->ry;
				RZ=atom->rz;
				WX=atom->wx;
				WY=atom->wy;
				WZ=atom->wz;
			}

			//Update the time
			t_delta += t_min;

			//Translate the particle to the moment it hits the BC
			stream_MD( atom,t_min );

			//Now the BC and the particle are touching. Let them collide.
			shiftBC_MD( shift,&WALL[chosenBC],atom );
			rotateBC_MD( &WALL[chosenBC],atom );
			//Find normal
			normal_MD( n,WALL[chosenBC],atom,DIM );
			//Normalize the normal vector
			norm( n,DIM );
			//Apply the BC to velocity
			velBC_MD( atom,&WALL[chosenBC],n,KBT );
			//Apply the BC to position
			posBC_MD( atom,WALL[chosenBC],n );
			rotatebackBC_MD( &WALL[chosenBC],atom );
			shiftbackBC( shift,&WALL[chosenBC] );
			//Update the time to stream for
			time = t_step - t_delta;
			//Let the BC stream the rest of the way
			stream_MD( atom,time );
			//Return to the top and try to move again for the rest of the time.
			W1 = calcW_MD(WALL[chosenBC],atom);
			if (W1 < -TOL){
				atom->rx=RX;
				atom->ry=RY;
				atom->rz=RZ;
				atom->wx=WX;
				atom->wy=WY;
				atom->wz=WZ;
				atom->vx=0.0;
				atom->vy=0.0;
				atom->vz=0.0;
			}
		}
	}
}

///
/// @brief	Checks if the MD particle is inside all the boundaries, if not it reports which 
/// 		boundary is crossed.
///
/// It checks the particle's position with respect to the BCs. If it is outside any of them 
/// it calculates the crosstime through crosstime_MD(). If the crosstime does not match the  
/// streaming timestep it raises a warning, but later in MD_BCcollision() routine the issue 
/// is solved.
///
/// @param WALL		All of the walls (BCs) that the particle might interact with.
/// @param atom		The MD particle.
/// @param t_min	The rest of the `time` remains for particle to stream after reducing the collision 
/// 				time.
/// @param chosenW 	It is used to determine if boundary conditions should be applied to the MD particle.
/// @param chosenBC The wall out of which the MD particle is.
/// @param time 	The total remaining time that the particle has in order to move.
/// @param t_step 	The MD timestep increment.
/// @see			chooseBC()
///
void chooseBC_MD( bc WALL[],particleMD *atom,double *t_min,double *chosenW,int *chosenBC,double time,double t_step ) {
	int i,flag;
	double t1,t2,tc;
	double tempW,shift[DIM];

	*t_min = time;
	*chosenW=1.;
	tempW = 1.;
	flag=0;

	for( i=0; i<NBC; i++ ) {
		//Shift BC due to periodic BCs
		shiftBC_MD( shift,&WALL[i],atom );
		rotateBC_MD( &WALL[i],atom );
		tempW=calcW_MD( WALL[i],atom );

		if( tempW < 0. ) {
			//Rewind the particle back to it's old position
			rewind_MD( atom,time );

			//Calculate crosstime
			crosstime_MD( atom,WALL[i],&t1,&t2,time );
			tc = chooseT( t_step,t1,t2,0,&flag );
			if( flag ) {
				printf( "Warning: Cross time unacceptable MD: %lf.\n",tc );
				// exit( 1 );		no need to exit. it is solved by further actions in MD_BCcollision
			}

			if( tc < *t_min ) {
				*t_min = tc;
				*chosenBC = i;
				*chosenW = tempW;
			}
			//Undo the rewind that's in this IF-statement
			stream_MD( atom,time );
		}

		//Shift BC back
		rotatebackBC_MD( &WALL[i],atom );
		shiftbackBC( shift,&WALL[i] );
	}
}

///
/// @brief	Determines if the BC must be shifted due to the periodicity of the control volume.
///   
/// It checks if the BC is periodic, then it calculates the shift and shifts the BC.
/// 
/// @param shift	This is how much the boundary must be shifted, gets calculated inside the routine.
/// @param WALL		One of the walls of the BCs that particle is interacting with.
/// @param atom 	The MD particle.
/// @see			shiftBC()
///
void shiftBC_MD( double *shift,bc *WALL,particleMD *atom ) {
	int k;

	for( k=0; k<_3D; k++ ) shift[k] =0.;
	//Don't shift planar surfaces
	if( WALL->P[0]>1.0 && WALL->P[1]>1.0 && WALL->P[2]>1.0 ) {
		//Check x-component
		if( atom->rx - WALL->Q[0] > 0.5*XYZ[0] ) shift[0] = XYZ[0];
		else if( atom->rx - WALL->Q[0] < -0.5*XYZ[0] ) shift[0] = -1.*XYZ[0];
		//Check y-component
		if( atom->ry - WALL->Q[1] > 0.5*XYZ[1] ) shift[1] = XYZ[1];
		else if( atom->ry - WALL->Q[1] < -0.5*XYZ[1] ) shift[1] = -1.*XYZ[1];
		//Check z-component
		if( atom->rz - WALL->Q[2] > 0.5*XYZ[2] ) shift[2] = XYZ[2];
		else if( atom->rz - WALL->Q[2] < -0.5*XYZ[2] ) shift[2] = -1.*XYZ[2];
	}
	//Shift BCs
	for( k=0; k<_3D; k++ ) WALL->Q[k] += shift[k];
}

///
/// @brief	Rotates the BC, if it has some orientation.
///
/// It rotates the particle's position, velocity, orientation about the BC surface instead using the
/// <a href="https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula">Rodrigues' rotation</a>. 
/// This routine does a forward rotation or backward rotation specified by a sign passed to it.
///
/// @param WALL		One of the walls of the BCs that the particle is interacting with.
/// @param atom		The MD particles.
/// @param sign 	The sign by which the orientation will happen.
/// @see 			rotateBC()
/// @see			rotatebackBC()
/// @see 			MPC_BCrotation()
///
void MD_BCrotation( bc *WALL,particleMD *atom, double sign ) {
	int i;
	double rotM[_3D][_3D];		//The rotation matrix
	double ax[_3D] = {1.0,0.0,0.0};	//x-axis
	double ay[_3D] = {0.0,1.0,0.0};	//y-axis
	double az[_3D] = {0.0,0.0,1.0};	//z-axis
	double oldQ[_3D] = {0.0,0.0,0.0};				//Particle position relative to the centre of the BC
	double newQ[_3D] = {0.0,0.0,0.0};				//Particle position relative to the centre of the BC
	double Q[_3D] = {atom->rx,atom->ry,atom->rz};	//Position of the atom
	double V[_3D] = {atom->vx,atom->vy,atom->vz};	//Velocity of the atom

	for( i=0; i<_3D; i++ ) oldQ[i] = Q[i]-WALL->Q[i];
	setRotMatrix3D( rotM,sign*WALL->O[0],sign*WALL->O[1],sign*WALL->O[2] );
	//Rotate the position into place
	dotprodMatVec( rotM,oldQ,newQ,DIM );
	atom->rx = WALL->Q[0] + newQ[0];
	atom->ry = WALL->Q[1] + newQ[1];
	atom->rz = WALL->Q[2] + newQ[2];
	//Rotate the velocity vector about each axis
	rodriguesRotation( V,az,sign*WALL->O[2] );
	rodriguesRotation( V,ay,sign*WALL->O[1] );
	rodriguesRotation( V,ax,sign*WALL->O[0] );
	atom->vx=V[0];
	atom->vy=V[1];
	atom->vz=V[2];
}

///
/// @brief	Checks if the BC has some orientation, if so it must be rotated.
///
/// MD_BCrotation() is used to do this in which the particle's position, velocity, orientation are
/// rotated about the BC surface instead. Uses NEGATIVE the angles since the particle is being rotated
/// instead of the BC.  
/// 
/// 
/// @param WALL		One of the walls of the BCs that the particle is interacting with.
/// @param atom 	The MD particle.
///	@note			The current implementation is very wasteful. Every \b particle
///					is rotated about the centre of each BC. While this is simplest, there are very many 
///					particles.
/// @see 			rotateBC()
///
void rotateBC_MD( bc *WALL,particleMD *atom ) {
	if(WALL->REORIENT) MD_BCrotation( WALL,atom,-1.0 );
}

///
/// @brief Undoes the rotateBC_MD().
/// 
/// @param WALL		One of the walls of the BCs that the particle is interacting with.
/// @param atom		The MD particle.
/// @see			rotatebackBC()
/// 
void rotatebackBC_MD( bc *WALL,particleMD *atom ) {
	if(WALL->REORIENT) MD_BCrotation( WALL,atom,1.0 );
}

///
/// @brief	Checks if the Boundary condition should be applied to the MD particle. 
///	
/// Calculates the distance of the particle from the wall, and based on the shape of the 
/// control volume it determines if the particle is outside/inside/on the wall.
///
/// @param WALL 	One of the walls of the BCs that the particle is interacting with.
/// @param atom 	The MD particle.
/// @return			The variable that detemines if the particle is inside(+)/outside(-)/on(0) the 
/// 				specified wall of the control volume.
/// @see			calcW()
///
double calcW_MD( bc WALL,particleMD *atom ){
	double terms, W=0.0;
	double Q[_3D];
	int i;

	Q[0] = atom->rx;
	Q[1] = atom->ry;
	Q[2] = atom->rz;

	if( feq(WALL.ROTSYMM[0],4.0) && feq(WALL.ROTSYMM[1],4.0) ) {
		for( i=0; i<_3D; i++ ) {
			terms = WALL.A[i] * ( Q[i]-WALL.Q[i] );
			if( WALL.ABS ) terms=fabs(terms);
			terms = smrtPow( terms,WALL.P[i] );
			W += terms;
		}
		terms = smrtPow( WALL.R,WALL.P[3] );
		if( WALL.ABS ) terms=fabs(terms);
		W -= terms;
	}
	else {
		printf( "Error:\tNon 4-fold symmetry not yet programmed for MD-BC interaction\n" );
		exit( 1 );
	}
	if( WALL.INV ) W *= -1.;

	return W;
}

///
/// @brief	The streaming step of the algorithm.
///
/// Using trans() routine translates the MD particle's position,
/// \f$ Q_{\mbox{New}} = Q_{\mbox{Old}} + t \times V \f$ 
/// in which \f$ Q_{\mbox{New}} \f$ is the new position of the particle,
/// \f$ Q_{\mbox{Old}} \f$ is the old position of the particle,
/// \f$ t \f$ is the streaming time and \f$ V \f$ is the velocity of the particle.
///
/// @param atom		The MD particle.
/// @param t		The time for which the particle must stream.
/// @note			No acceleration during time `t`.
/// @see			stream_P()
/// @see            trans()
///
void stream_MD( particleMD *atom,double t ) {
	atom->rx = trans( t,atom->vx,atom->rx );
	atom->ry = trans( t,atom->vy,atom->ry );
	atom->rz = trans( t,atom->vz,atom->rz );
}

///
/// @brief	Rewinds the particle back to its old position.
///
/// Using rewind_trans() brings back the particle to its old position,
/// \f$ Q_{\mbox{New}} = Q_{\mbox{Old}} - t \times V \f$ in which \f$ Q_{\mbox{New}} \f$
/// is the new position of the particle,
/// \f$ Q_{\mbox{Old}} \f$ is the old position of the particle, \f$ t \f$ is the streaming time
/// and \f$ V \f$ is the velocity of the particle.
///
/// @param atom		The MD particle.
/// @param time		The time for which the particle streams backward.
/// @see			rewind_P()
/// @see			rewind_trans()
///
void rewind_MD( particleMD *atom,double time ) {
	atom->rx = rewind_trans(time,atom->vx,atom->rx);
	atom->ry = rewind_trans(time,atom->vy,atom->ry);
	atom->rz = rewind_trans(time,atom->vz,atom->rz);
}

///
/// @brief	Calculates when the MD particle crosses the BC.
///
/// It calculates the time takes for the particle to cross the boundary by solving the trajectory 
///	equation, \f$ \left[ \left(Q - Q_c \right) + V \times t \right]^2 = R ^2 \f$, in which \f$ Q \f$
/// is the position of the particle, \f$ Q_c \f$ is the position of the center of the control volume,\f$ 
/// V \f$ is the velocity of the particle,
/// \f$ t \f$ is the streaming time, and \f$ R \f$ is the radius of the control volume.
/// It can also use the <a href="https://en.wikipedia.org/wiki/Secant_method">secant method</a>.
///  
/// @param atom		The MD particle.
/// @param WALL		One of the walls of the BCs that the particle is interacting with.
/// @param tc_pos	One of the crosstimes.
/// @param tc_neg	One of the crosstimes.
/// @param t_step	The maximum streaming time.
///	@see			secant_time_MD(), crosstime() 
/// @note			In case the wall is \b planar or the \b the \b secant \b method is used
/// 				\f$ tc_pos = tc_neg \f$.
///
void crosstime_MD( particleMD *atom,bc WALL,double *tc_pos, double *tc_neg,double t_step ) {
	double a=0.0,b=0.0,c=0.0;

	// Planar Wall
	if( WALL.PLANAR || ( feq(WALL.P[0],1.0) && feq(WALL.P[1],1.0) && feq(WALL.P[2],1.0) && feq(WALL.P[3],1.0) )) {
		*tc_pos = WALL.R;
		//x-component
		*tc_pos += WALL.A[0]*(WALL.Q[0]-atom->rx);
		//y-component
		*tc_pos += WALL.A[1]*(WALL.Q[1]-atom->ry);
		//z-component
		*tc_pos += WALL.A[2]*(WALL.Q[2]-atom->rz);

		*tc_pos /= (WALL.A[0]*atom->vx + WALL.A[1]*atom->vy + WALL.A[2]*atom->vz);
		//There is only one time.
		*tc_neg = *tc_pos;
	}
	// Ellipsoid
	else if((DIM == 2 && feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0)) || (DIM > 2 && feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0) && feq(WALL.P[2],2.0))) {
		//x-component
		a += WALL.A[0]*WALL.A[0]*atom->vx*atom->vx;
		b += WALL.A[0]*WALL.A[0]*atom->vx*(atom->rx-WALL.Q[0]);
		c += WALL.A[0]*WALL.A[0]*(atom->rx*atom->rx-2.*atom->rx*WALL.Q[0]+WALL.Q[0]*WALL.Q[0]);
		//y-component
		a += WALL.A[1]*WALL.A[1]*atom->vy*atom->vy;
		b += WALL.A[1]*WALL.A[1]*atom->vy*(atom->ry-WALL.Q[1]);
		c += WALL.A[1]*WALL.A[1]*(atom->ry*atom->ry-2.*atom->ry*WALL.Q[1]+WALL.Q[1]*WALL.Q[1]);
		//z-component
		a += WALL.A[2]*WALL.A[2]*atom->vz*atom->vz;
		b += WALL.A[2]*WALL.A[2]*atom->vz*(atom->rz-WALL.Q[2]);
		c += WALL.A[2]*WALL.A[2]*(atom->rz*atom->rz-2.*atom->rz*WALL.Q[2]+WALL.Q[2]*WALL.Q[2]);

		b *= 2.0;
		c -= smrtPow(WALL.R,WALL.P[3]);
		//Use the quadratic formula
		*tc_neg = - b - sqrt(b*b-4.*a*c);
		*tc_neg /= (2.*a);
		*tc_pos = - b + sqrt(b*b-4.*a*c);
		*tc_pos /= (2.*a);
	}
	else {
		//Must use secant method to determine cross times
		*tc_pos = secant_time_MD( atom,WALL,t_step );
		*tc_neg = *tc_pos;
	}
}

///
/// @brief	Numerically determines the crossing times.
///
/// It uses the <a href="https://en.wikipedia.org/wiki/Secant_method">secant method</a> to calculate the 
/// the cross time. The secant method is a root-finding algorithm that uses a succession of roots of 
/// secant lines to better approximate a root of a function f, which in this case is the trajectory of
/// the particle.
///
/// @param atom		The MD particle.
/// @param WALL		One of the walls of the BCs that the particle is interacting with.
/// @param t_step	The maximum streaming time.
/// @return 		The crosstime.
///
double secant_time_MD( particleMD *atom,bc WALL,double t_step ) {
	double Qi[DIM],QiM1[DIM];
	double ti,tiM1,root;
	double fi,fiM1;
	int iter=0;

	//Rewind the particle back to it's old position
	//x-component
	QiM1[0] = atom->rx;
	Qi[0] = trans(t_step,atom->vx,atom->rx);
	//1-component
	QiM1[1] = atom->ry;
	Qi[1] = trans(t_step,atom->vy,atom->ry);
	//z-component
	QiM1[2] = atom->rz;
	Qi[2] = trans(t_step,atom->vz,atom->rz);

	ti = 0.;
	tiM1 = t_step;
	root = ti;

	//Secant Loop
	do {
		iter++;
		// Calculate the surface function for the particles' positions at these times
		fi = surf_func( WALL,Qi,DIM );
		fiM1 = surf_func( WALL,QiM1,DIM );

		root = ti - fi * ( ti - tiM1 )/ ( fi - fiM1 );
		tiM1 = ti;
		ti = root;
		// Calculate the particles' positions at these times
		//x-component
		QiM1[0] = trans(tiM1,atom->vx,atom->rx);
		Qi[0] = trans(ti,atom->vx,atom->rx);
		//y-component
		QiM1[1] = trans(tiM1,atom->vy,atom->ry);
		Qi[1] = trans(ti,atom->vy,atom->ry);
		//z-component
		QiM1[2] = trans(tiM1,atom->vz,atom->rz);
		Qi[2] = trans(ti,atom->vz,atom->rz);

	} while( fabs(ti-tiM1) > TOL && (fabs(fi-fiM1) > TOL) );
	return root;
}

///
/// @brief	Finds the normal to the surface that the particle is
///			presently ON.
///
/// It takes the gradient of \f$ ( a(x-h) )^p + (b(y-k))^p + (c(z-l))^p - r = 0 \f$, that is the equation
/// for the control volume, since the gradient is equal to the normal. For powers of 1 and 2 it takes the
/// shortcuts and uses the specific solution programmed in. For higher powers it uses a more general 
/// solution. 
/// 
/// @param	n			The normal vector to the surface.
/// @param WALL 		One of the walls of the BCs that the particle is interacting with.
/// @param atom 		The MD particle.
/// @param dimension	The dimenson of the control volume.
/// @return				The normal vector to the surface.
/// @see 				normal()
///
double *normal_MD( double *n,bc WALL,particleMD *atom,int dimension ) {
	int i;

	if( WALL.PLANAR || ( feq(WALL.P[0],1.0) && feq(WALL.P[1],1.0) && feq(WALL.P[2],1.0) && feq(WALL.P[3],1.0))) {
		for( i=0; i<dimension; i++ ) n[i] = WALL.A[i];
	}
	// Issue: for 2D, can't keep feq(WALL.P[2], 2.0) as P[2]=1.0.
	else if( feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0) && feq(WALL.P[2],2.0) ) {
		n[0] = 2.*WALL.A[0]*(atom->rx-WALL.Q[0]);
		n[1] = 2.*WALL.A[1]*(atom->ry-WALL.Q[1]);
		n[2] = 2.*WALL.A[2]*(atom->rz-WALL.Q[2]);
	}
	else {
		n[0] = (WALL.P[0]) *smrtPow( WALL.A[0]*(atom->rx-WALL.Q[0]) , WALL.P[0]-1.0 );
		n[1] = (WALL.P[1]) *smrtPow( WALL.A[1]*(atom->ry-WALL.Q[1]) , WALL.P[1]-1.0 );
		n[2] = (WALL.P[2]) *smrtPow( WALL.A[2]*(atom->rz-WALL.Q[2]) , WALL.P[2]-1.0 );
	}

	return n;
}

///
/// @brief	This subroutine applies the BC transformation to the velocity of the MD particle.
///
///	It transforms the velocity (the normal and tangential components) of the MD particle considering the 
/// conditions at the surface of the boundary. For instance, it can \b conserve the \b energy/momentum/
/// \b angular \b momentum using impulse method or it can apply the rule method such as \b bounceback or 
/// \b reflection or \b periodic which does NOT necesarily conserve momentum. The BCs global variables 
/// defined in the definition.h set how particle velosity will be updated.  
///
/// @param atom		The MD particle.
/// @param WALL 	One of the walls of the BCs that the particle is interacting with.
/// @param n 		The normal vector to the surface of the wall.
/// @param KBT 		Thermal energy.
/// @see			velBC()
///
void velBC_MD( particleMD *atom,bc *WALL,double n[_3D],double KBT ) {
	double V[_3D],VN[_3D],VT[_3D],VR[_3D],R[_3D],zip[_3D];
	double atom_POS[_3D],atom_VEL[_3D];
	double IIpart[_3D][_3D],IIwall[_3D][_3D];
	double IMpart,IMwall;		//Inverse mass
	double J=1.;			//Impulse
	double rand;			//Random number
	int i,j;

	//Zero everything
	for( i=0; i<_3D; i++ ) {
		V[i] = 0.;
		VN[i] = 0.;
		VT[i] = 0.;
		VR[i] = 0.;
		R[i] = 0.;
		zip[i] = 0.;
	}

	atom_POS[0] = atom->rx;
	atom_POS[1] = atom->ry;
	atom_POS[2] = atom->rz;
	atom_VEL[0] = atom->vx;
	atom_VEL[1] = atom->vy;
	atom_VEL[2] = atom->vz;

	//Calculate rotational velocity of BC
	for( i=0; i<DIM; i++ ) R[i] = atom_POS[i] - WALL->Q[i];
	crossprod( WALL->L,R,VR );
	//Begin Determining the direction of the impulse
	for( i=0; i<DIM; i++ ) V[i] = atom_VEL[i] - WALL->V[i] - VR[i];
	//Calculate normal and tangential components of velocity
	proj( V,n,VN,DIM );
	tang( V,VN,VT,DIM );

	//Point particles do not have angular momentum
	for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) IIpart[i][j] = 0.;
	IMpart = 1. / atom->mass;
	if( !(WALL->DSPLC)) {
		for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) IIwall[i][j] = 0.;
		IMwall = 0.;
	}
	else if ( WALL->DSPLC ) {
		invert3x3(IIwall,WALL->I);
		IMwall = 1.0 / WALL->MASS;
	}
	else {
		printf( "Error: WALL.DSPLC must be 0 or 1.\n" );
		exit( 1 );
	}

	//The energy/momentum/angular momentum conserving impulse method
	if( WALL->COLL_TYPE == BC_IMP ) {
		//Transform the normal and tangential components
		for( i=0; i<DIM; i++ ) {
			VN[i] *= WALL->MVN;
			VN[i] += WALL->DVN*n[i];
			VT[i] *= WALL->MVT;
			VT[i] += WALL->DVT;
		}
		//Combine normal and tangential components
		V[0] = VN[0] + VT[0];
		V[1] = VN[1] + VT[1];
		if( DIM >= _3D ) V[2] = VN[2] + VT[2];
		//The impulse's direction is the difference between the two (final minus initial)
		for( i=0; i<DIM; i++ ) V[i] -= ( atom_VEL[i] - WALL->V[i] - VR[i] );
		//Normalize V. We only want the unit vector. Conservation will determine its magnitude
		norm( V,DIM );

		J = impulse( V,atom_VEL,WALL->V,atom_POS,WALL->Q,zip,WALL->L,IMpart,IMwall,IIpart,IIwall,atom_POS,WALL->E );
	}
	//The rule method such as bounceback or reflection or periodic which does NOT necesarily conserve momentum
	else if( WALL->COLL_TYPE == BC_SURF_RULES ) {
		//Transform the normal and tangential components
		for( i=0; i<DIM; i++ ) {
			VN[i] *= WALL->MVN;
			VN[i] += WALL->DVN*n[i];
			VT[i] *= WALL->MVT;
			VT[i] += WALL->DVT;
		}
		//Combine normal and tangential components
		V[0] = VN[0] + VT[0];
		V[1] = VN[1] + VT[1];
		if( DIM >= _3D ) V[2] = VN[2] + VT[2];
		//Move the velocity out of the particle's rest frame and back into the lab frame
		for( i=0; i<DIM; i++ ) V[i] += WALL->V[i] + VR[i];
		//Make V the change in momentum instead of the velocity by subtracting the old velocity and times mass
		for( i=0; i<DIM; i++ ) {
			V[i] -= atom_VEL[i];
			V[i] *= atom->mass;
		}
		//Set J=1 because the total change in momentum (i.e. impulse) is in V.
		J = 1.;
	}
	else if( WALL->COLL_TYPE == BC_THERMO_SURF ) {
		//Transform the normal and tangential components
		for( i=0; i<DIM; i++ ) {
			VN[i] *= WALL->MVN;
			VN[i] += WALL->DVN*n[i];
			VT[i] *= WALL->MVT;
			VT[i] += WALL->DVT;
		}
		norm( VN,DIM );
		norm( VT,DIM );
		// Apply Thermo boundary conditions
		rand = genrand_gaussMB( WALL->KBT,atom->mass );
		for( i=0; i<_3D; i++ ) VT[i] *= rand;
		rand = genrand_rayleigh( sqrt(WALL->KBT/atom->mass) );
		for( i=0; i<_3D; i++ ) VN[i] *= rand;
		//Combine normal and tangential components
		V[0] = VN[0] + VT[0];
		V[1] = VN[1] + VT[1];
		if( DIM >= _3D ) V[2] = VN[2] + VT[2];

		/************************************************************************/
		/************************************************************************/
		//Needs a thermostat to be on!
		/************************************************************************/
		/************************************************************************/
		//Move the velocity out of the particle's rest frame and back into the lab frame
		for( i=0; i<DIM; i++ ) V[i] += WALL->V[i] + VR[i];
		//Make V the change in momentum instead of the velocity by subtracting the old velocity and times mass
		for( i=0; i<DIM; i++ ) {
			V[i] -= atom_VEL[i];
			V[i] *= atom->mass;
		}
		//Set J=1 because the total change in momentum (i.e. impulse) is in V.
		J = 1.;
	}
	else if( WALL->COLL_TYPE == BC_HALF_RULES ) {
		printf("Error: Rule-based solvent/BC collisions at t/2 not yet operational.\n");
		exit(1);
	}
	else if( WALL->COLL_TYPE == BC_THERMO_HALF ){
		printf("Error: Probabilistic solvent/BC collisions at t/2 not yet operational.\n");
		exit(1);
	}
	else {
		printf("Error: Solvent/BC Collision type unrecognized.\n");
		exit(1);
	}

	//Use the impulse to set the velocity of particleMPC
	atom->vx += V[0] * J * IMpart;
	atom->vy += V[1] * J * IMpart;
	atom->vz += V[2] * J * IMpart;
	if( WALL->DSPLC ) {
		//Set the velocity of BC
		for( i=0; i<DIM; i++ ) WALL->V[i] -= V[i] * J * IMwall;
		//Set the angular velocity of BC
		//Since VN isn't being used, use VN as the difference
		for( i=0; i<_3D; i++ ) VN[i] = atom_POS[i] - WALL->Q[i];
		//Since VT isn't being used, use VT as the crossprod.
		crossprod( VN,V,VT );
// 		dotprodmat( VT,IIwall,VN,_3D );
		dotprodMatVec( IIwall,VT,VN,_3D );
		for( i=0; i<_3D; i++) WALL->dL[i] -= VN[i] * J;
	}
}

///
/// @brief	This subroutine applies the BC transformation to the position of the MD particle. 
///
/// It updates the position of the MD particle by applying the normal and tangential displacements 
/// specified in the input file, when crossing a periodic boundary.
///  
/// @param atom		The MD particle.
/// @param WALL 	One of the walls of the BCs that the particle is interacting with.
/// @param n 		The normal vector to the surface of the wall.
/// @note			It does NOT stream! That is done in a seperate routine.
/// @see			posBC()
///	
void posBC_MD( particleMD *atom,bc WALL,double n[_3D] ) {
	double PN[_3D],PT[_3D],temp[_3D];
	int i;

	temp[0] = atom->rx;
	temp[1] = atom->ry;
	temp[2] = atom->rz;

	proj( temp,n,PN,DIM );		//Calculate normal component of the position
	tang( temp,PN,PT,DIM );		//Calculate tangential component of position

	//Transform the position
	for( i=0; i<DIM; i++ ) {
		PN[i] += WALL.DN * n[i];//Transform normal component
		PT[i] += WALL.DT;		//Transform tangential component
	}
	//Combine normal and tangential components
	atom->rx = PN[0] + PT[0];
	if( DIM > _1D ){
		atom->ry = PN[1] + PT[1];
	}
	if( DIM > _2D ){
		atom->rz = PN[2] + PT[2];
	}
	
}

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ********** SWIMMERS PASSING BCs ********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */

///
/// @brief	Applies the BCs to the MD particle forming the swimmer body, in case it is required.
///
/// It checks if the particle is inside the boundaries. If it's inside the boundaries it does not
/// do anything. But if it is not, it rewinds it to its old position, then it calculates the time taken
/// for the particle to collide with the boundary all through chooseBC_swimmer(). it streams that time, 
/// collides with the boundary, then it streams in to the control volume for the rest of the streaming 
/// time. If the function fails to bring it inside the control volume, it is simply brought back to 
/// its old position with no velocity.
///
/// @param atom		The MD particle, being either the head or the middle monomer of the swimmer.
/// @param WALL 	One of the walls of the BCs that the particle is interacting with.
/// @param SS 		It specifies the type (features) of the swimmer to which the monomer belongs.
/// @param t_step	The MD timestep increment.
/// @see			MPC_BCcollision()
///  
void swimmer_BCcollision( smono *atom,bc WALL[],specSwimmer SS,double t_step ) {
	double t_delta;			//time passed so far
	double time;			//time left to move for
	double t_min=0.;		//smallest time
	int chosenBC=0;			//Particle to go with t_min
	int flag=1;			//flag for if should keep looping. Loop while 1.
	int cnt=0;				//Count the loop iterations
	double n[_3D] = {0.,0.,0.};	//Normal to the surface
	double W,W1=0.0;
	double shift[DIM];
	double RX=0.0,RY=0.0,RZ=0.0;

	t_delta = 0.;
	time = t_step;
	while( flag ) {
		//We must check if the particle is inside any of the BCs
		chooseBC_swimmer( WALL,atom,&t_min,&W,&chosenBC,time,t_step );

		//If no particles were inside then we are done.
		if( W > TOL ) flag = 0;
		//Otherwise, COLLISON
		else {
			cnt++;
			// #ifdef DBG
			// 	if( DBUG == DBGSWIMMER ) printf( "\tW=%f BC=%d\n",W,chosenBC );
			// #endif
			//We have the BC to collide with and the time at which it collided
			//Rewind the particle back to it's old position
			rewind_swimmer( atom,time );
			if(cnt==1){
				RX=atom->Q[x_];
				RY=atom->Q[y_];
				RZ=atom->Q[z_];
			}
			//Update the time
			t_delta += t_min;

			//Translate the particle to the moment it hits the BC
			stream_swimmer( atom,t_min );

			//Now the BC and the particle are touching. Let them collide.
			shiftBC_swimmer( shift,&WALL[chosenBC],atom );
			rotateBC_swimmer( &WALL[chosenBC],atom );
			//Find normal
			normal_swimmer( n,WALL[chosenBC],atom,DIM );
			//Normalize the normal vector
			norm( n,DIM );
			//Apply the BC to velocity
			velBC_swimmer( atom,&WALL[chosenBC],SS,n );
			//Apply the BC to position
			posBC_swimmer( atom,WALL[chosenBC],n );
			rotatebackBC_swimmer( &WALL[chosenBC],atom );
			shiftbackBC( shift,&WALL[chosenBC] );

			//Update the time to stream for
			time = t_step - t_delta;
			//Let the BC stream the rest of the way
			stream_swimmer( atom,time );
			//Return to the top and try to move again for the rest of the time.
			W1 = calcW_swimmer(WALL[chosenBC],atom);
			if (W1 < 0.0){
				atom->Q[x_]=RX;
				atom->Q[y_]=RY;
				atom->Q[z_]=RZ;
				atom->V[x_]=0.0;
				atom->V[y_]=0.0;
				atom->V[z_]=0.0;
			}
		}
	}
}

///
/// @brief	Checks if the MD particle forming the swimmer is inside all the boundaries, if not
///			it reports the boundary.
///
/// It checks the particle's position with respect to the BCs. If it is outside any of them 
/// it calculates the crosstime through crosstime_swimmer(). If the crosstime does not match the
/// streaming timestep it raises a warning, but later in swimmer_BCcollision routine the issue is solved.
///
/// @param WALL 	All of the walls (BCs) that the particle might interact with.
/// @param atom 	The MD particle, being either the head or the middle monomer of the swimmer.
/// @param t_min 	The rest of the `time` remains for the particle to stream after reducing the 		collision time.
/// @param chosenW 	It is used to determine if the boundary conditions should be applied to the MD particle.
/// @param chosenBC The wall out of which the MD particle is.
/// @param time 	The total remaining time that the particle has to move.
/// @param t_step 	The MD timestep increment.
/// @see			chooseBC()
///
void chooseBC_swimmer( bc WALL[],smono *atom,double *t_min,double *chosenW,int *chosenBC,double time,double t_step ) {
	int i,flag;
	double t1,t2,tc;
	double tempW,shift[DIM];

	*t_min = time;
	*chosenW=1.;
	tempW = 1.;
	flag=0;

	for( i=0; i<NBC; i++ ) {
		//Shift BC due to periodic BCs
		shiftBC_swimmer( shift,&WALL[i],atom );
		rotateBC_swimmer( &WALL[i],atom );
		tempW=calcW_swimmer( WALL[i],atom );

		if( tempW < 0. ) {
			//Rewind the particle back to it's old position
			rewind_swimmer( atom,time );

			//Calculate crosstime
			crosstime_swimmer( atom,WALL[i],&t1,&t2,time );
			tc = chooseT( t_step,t1,t2,0,&flag );
			if( flag ) {
				printf( "Warning: Cross time unacceptable swimmers: %lf.\n",tc );
				// exit( 1 );      //no need to exit, it is solved by further actions in swimmer_BCcollision
			}

			if( tc < *t_min ) {
				*t_min = tc;
				*chosenBC = i;
				*chosenW = tempW;
			}
			//Undo the rewind that's in this IF-statement
			stream_swimmer( atom,time );
		}

		//Shift BC back
		rotatebackBC_swimmer( &WALL[i],atom );
		shiftbackBC( shift,&WALL[i] );
	}
}

///
/// @brief	Determines if the BC must be shifted due to the periodicity of the control volume.
///
/// It checks if the BC is periodic, then it calculates the shift and shifts the BC.
///
/// @param shift	This is how much the boundary must be shifted, gets calculated inside the routine.
/// @param WALL 	One of the walls of the BCs that the particle is interacting with.
/// @param atom 	The MD particle, being either the head or the middle monomer of the swimmer. 
///
void shiftBC_swimmer( double *shift,bc *WALL,smono *atom ) {
	int d;

	for( d=0; d<_3D; d++ ) shift[d] =0.0;
	//Don't shift planar surfaces
	if( WALL->P[0]>1.0 && WALL->P[1]>1.0 && WALL->P[2]>1.0 ) for( d=0; d<DIM; d++ ) {
		if( atom->Q[d] - WALL->Q[d] > 0.5*XYZ[d] ) shift[d] = XYZ[d];
		else if( atom->Q[d] - WALL->Q[d] < -0.5*XYZ[d] ) shift[d] = -1.0*XYZ[d];
	}
	//Shift BCs
	for( d=0; d<DIM; d++ ) WALL->Q[d] += shift[d];
}

///
/// @brief	Rotates the BC, if it has some orientation.
///
/// It rotates the particle's posotion, velocity, orientation about the BC surface instead using the 
/// <a href="https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula">Rodrigues' rotation</a>.
/// This routine does the forward rotation or backward rotation by specifying a sign passed to it.
///
/// @param WALL		One of the walls of the BCs that the particle is interacting with.
/// @param atom 	The MD particles, being either the head or the middle monomer of the swimmer.
/// @param sign 	The sign by which the orientation will happen.
/// @see			rotateBC()
/// @see			rotatebackBC()
/// @see 			MPC_BCrotation()
///
void swimmer_BCrotation( bc *WALL,smono *atom, double sign ) {
	int i;
	double rotM[_3D][_3D];		//The rotation matrix
	double ax[_3D] = {1.0,0.0,0.0};	//x-axis
	double ay[_3D] = {0.0,1.0,0.0};	//y-axis
	double az[_3D] = {0.0,0.0,1.0};	//z-axis
	double oldQ[_3D] = {0.0,0.0,0.0};				//Particle position relative to the centre of the BC
	double newQ[_3D] = {0.0,0.0,0.0};				//Particle position relative to the centre of the BC

	for( i=0; i<DIM; i++ ) oldQ[i] = atom->Q[i]-WALL->Q[i];
	if( DIM>_2D ) setRotMatrix3D( rotM,sign*WALL->O[0],sign*WALL->O[1],sign*WALL->O[2] );
	else setRotMatrix2D( rotM,sign*WALL->O[2] );
	//Rotate the position into place
	dotprodMatVec( rotM,oldQ,newQ,DIM );
	for( i=0; i<DIM; i++ ) atom->Q[i] = WALL->Q[i] + newQ[i];
	//Rotate the velocity vector about each axis
	//Z-axis
	rodriguesRotation( atom->V,az,sign*WALL->O[2] );
	//X- & Y-axes
	if( DIM>_2D ) {
		rodriguesRotation( atom->V,ay,sign*WALL->O[1] );
		rodriguesRotation( atom->V,ax,sign*WALL->O[0] );
	}
}

///
/// @brief	Checks if the BC has some orientation, if so it must be rotated.
///
/// swimmer_BCrotation() is used to do this in which the particle's position, velocity, orientation are 
/// rotated about the BC surface instead. Uses NEGATIVE the angles since the particle is being rotated
/// instead of the BC.
///
/// @param WALL		One of the walls of the BCs that the particle is interacting with.
/// @param atom 	The MD particle, being either the head or the middle monomer of the swimmer.
///	@note			The current implementation is very wasteful. Every \b particle is rotated
///					about the centre of each BC. While this is simplest, there are very many 
///					particles.
/// @see			rotateBC()
///
void rotateBC_swimmer( bc *WALL,smono *atom ) {
	if(WALL->REORIENT) swimmer_BCrotation( WALL,atom,-1.0 );
}

///
/// @brief 	Undoes the rotateBC_swimmer().
/// 
/// @param WALL		One of the walls of the BCs that the article is interacting with.
/// @param atom		The MD particle, being either the head or the middle monomer of the swimmer.
/// @see			rotatebackBC()
///
void rotatebackBC_swimmer( bc *WALL,smono *atom ) {
	if(WALL->REORIENT) swimmer_BCrotation( WALL,atom,1.0 );
}

///
/// @brief	Checks if the boundary condition should be applied to the MD particle forming the swimmer. 
///	
/// Calculates the distance of the particle from the wall, and based on the shape of the 
/// control volume determines if the particle is outside/inside/on the wall.
///
/// @param WALL 	One of the walls of the BCs that the particle is interacting with.
/// @param atom 	The MD particle, being either the head or the middle monomer of the swimmer.
/// @return			The variable that detemines if the particle is inside(+)/outside(-)/on(0) the 
///  				specified wall of the control volume.
/// @see			calcW()
///
double calcW_swimmer( bc WALL,smono *atom ) {
	double terms, W=0.0;
	int d;

	if( feq(WALL.ROTSYMM[0],4.0) && feq(WALL.ROTSYMM[1],4.0) ) {
		for( d=0; d<DIM; d++ ) {
			terms = WALL.A[d] * ( atom->Q[d]-WALL.Q[d] );
			if( WALL.ABS ) terms=fabs(terms);
			terms = smrtPow( terms,WALL.P[d] );
			W += terms;
		}
		terms = smrtPow( WALL.R,WALL.P[3] );
		if( WALL.ABS ) terms=fabs(terms);
		W -= terms;
	}
	else {
		printf( "Error:\tNon 4-fold symmetry not yet programmed for swimmer-BC interactions\n" );
		exit( 1 );
	}
	if( WALL.INV ) W *= -1.;
	return W;
}

///
/// @brief	The streaming step of the algorithm, swimmers' version.
///
/// Using trans() routine translates the MD particle's position, 
/// \f$ Q_{\mbox{New}} = Q_{\mbox{Old}} + t \times V \f$, in which \f$ Q_{\mbox{New}} \f$ is the new 
/// position of the particle,
/// \f$ Q_{\mbox{Old}} \f$ is the old position of the particle, \f$ t \f$ is the streaming time and \f$ 
/// V \f$ is the velocity of the particle.
///
/// @param atom		The MD particle, being either the head or the middle monomer of the swimmer.
/// @param t		The time for which the particle must stream.
/// @note			No acceleration during time `t`.
/// @see			stream_P()
/// @see			trans()
///
void stream_swimmer( smono *atom,double t ) {
    int d;
    for( d=0; d<DIM; d++ ) atom->Q[d] = trans( t,atom->V[d],atom->Q[d] );
}

///
/// @brief	Rewinds the particle back to its old position, swimmers' version.
///
/// Using rewind_trans() brings back the particle to its old position,
/// \f$ Q_{\mbox{New}} = Q_{\mbox{Old}} - t \times V \f$, in which \f$ Q_{\mbox{New}} \f$ is the new 
/// position of the particle, \f$ Q_{\mbox{Old}} \f$ is the old position of the particle,
/// \f$ t \f$ is the streaming time and \f$ V \f$ is the velocity of the particle.
/// 
/// @param atom		The MD particle, being either the head or the middle monomer of the swimmer.
/// @param time 	The time for which the particle streams backward.
/// @see			rewind_P()
/// @see			rewind_trans()
///
void rewind_swimmer( smono *atom,double time ) {
	int d;
    for( d=0; d<DIM; d++ ) atom->Q[d] = rewind_trans( time,atom->V[d],atom->Q[d] );
}

///
/// @brief	Calculates when the MD particle crosses the BC, swimmers' version.
///
/// It calculates the time taken for the particle to cross the boundary either by solving the trajectory 
///	equation, \f$ \left[ \left(Q - Q_c \right) + V \times t \right]^2 = R ^2 \f$, in which \f$ Q \f$
/// is the position of the particle, \f$ Q_c \f$ is the position of the center of the control volume,\f$ 
/// V \f$ is the velocity of the particle,
/// \f$ t \f$ is the streaming time, \f$ R \f$ is the radius of the control volume. 
/// It can also use the <a href="https://en.wikipedia.org/wiki/Secant_method">secant method</a>.
/// 
/// @param atom		The MD particle, being either the head or the middle monomer of the swimmer.
/// @param WALL 	One of the walls of the BCs that the particle is interacting with.
/// @param tc_pos 	One of the crosstimes.
/// @param tc_neg 	One of the crosstimes.
/// @param t_step 	The maximum streaming time.
///	@see			secant_time_swimmer()
/// @see 			crosstime()
/// @note			In case the wall is \b planar or the \b the \b secant \b method is used
/// 				\f$ tc_pos = tc_neg \f$.
///
void crosstime_swimmer( smono *atom,bc WALL,double *tc_pos, double *tc_neg,double t_step ) {
    int d;
	double a=0.0,b=0.0,c=0.0;

	// Planar Wall
	if( WALL.PLANAR || ( feq(WALL.P[0],1.0) && feq(WALL.P[1],1.0) && feq(WALL.P[2],1.0) && feq(WALL.P[3],1.0) )) {
		*tc_pos = WALL.R;
		for( d=0; d<DIM; d++ ) *tc_pos += WALL.A[d]*(WALL.Q[d]-atom->Q[d]);
		*tc_pos /= (WALL.A[0]*atom->V[0] + WALL.A[1]*atom->V[1] + WALL.A[2]*atom->V[2]);
		//There is only one time.
		*tc_neg = *tc_pos;
	}
	// Ellipsoid
	else if((DIM == 2 && feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0)) || (DIM > 2 && feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0) && feq(WALL.P[2],2.0))) {
		for( d=0; d<DIM; d++ ) {
			a += WALL.A[d]*WALL.A[d]*atom->V[d]*atom->V[d];
			b += WALL.A[d]*WALL.A[d]*atom->V[d]*(atom->Q[d]-WALL.Q[d]);
			c += WALL.A[d]*WALL.A[d]*(atom->Q[d]*atom->Q[d]-2.0*atom->Q[d]*WALL.Q[d]+WALL.Q[d]*WALL.Q[d]);
		}
		b *= 2.0;
		c -= smrtPow(WALL.R,WALL.P[3]);
		//Use the quadratic formula
		*tc_neg = - b - sqrt(b*b-4.*a*c);
		*tc_neg /= (2.*a);
		*tc_pos = - b + sqrt(b*b-4.*a*c);
		*tc_pos /= (2.*a);
	}
	else {
		//Must use secant method to determine cross times
		*tc_pos = secant_time_swimmer( atom,WALL,t_step );
		*tc_neg = *tc_pos;
	}
}

///
/// @brief	Numerically determines the crossing times, swimmers' version.
///
/// It uses the <a href="https://en.wikipedia.org/wiki/Secant_method">secant method</a> to calculate the 
/// the cross time. The secant method is a root-finding algorithm that uses a succession of roots of 
/// secant lines to better approximate a root of a function f, which in this case is the trajectory of
/// the particle.
/// 
/// @param atom		The MD particle, being either the head or the middle monomer of the swimmer.
/// @param WALL 	One of the walls of the BCs that the particle is interacting with.
/// @param t_step 	The maximum streaming time.
/// @return 		The crosstime.
///
double secant_time_swimmer( smono *atom,bc WALL,double t_step ) {
	double Qi[DIM],QiM1[DIM];
	double ti,tiM1,root;
	double fi,fiM1;
	int iter=0;
	int d;

	//Rewind the particle back to it's old position
	for( d=0; d<DIM; d++ ) {
		QiM1[d] = atom->Q[d];
		Qi[d] = trans(t_step,atom->V[d],atom->Q[d]);
	}
	ti = 0.0;
	tiM1 = t_step;
	root = ti;

	//Secant Loop
	do {
		iter++;
		// Calculate the surface function for the particles' positions at these times
		fi = surf_func( WALL,Qi,DIM );
		fiM1 = surf_func( WALL,QiM1,DIM );

		root = ti - fi * ( ti - tiM1 )/ ( fi - fiM1 );
		tiM1 = ti;
		ti = root;
		// Calculate the particles' positions at these times
		for( d=0; d<DIM; d++ ) {
			QiM1[d] = trans(tiM1,atom->V[d],atom->Q[d]);
			Qi[d] = trans(ti,atom->V[d],atom->Q[d]);
		}
	} while( fabs(ti-tiM1) > TOL && (fabs(fi-fiM1) > TOL) );
	return root;
}

///
/// @brief	Finds the normal to the surface that the particle is
///			presently ON.
///
/// It takes the gradient of \f$ ( a(x-h) )^p + (b(y-k))^p + (c(z-l))^p - r = 0 \f$, that is the equation
/// for the control volume, since the gradient is equal to the normal. For powers of 1 and 2 it takes  
/// the shortcuts and uses the specific solution programmed in. For higher powers it uses a more general 
/// solution. 
///
/// @param n 			The normal vector to the surface.
/// @param WALL 		One of the walls of the BCs that the particle is interacting with.
/// @param atom 		The MD particle, being either the head or the middle monomer of the swimmer.
/// @param dimension	The dimenson of the control volume.
/// @return				The normal vector to the surface.
/// @see				normal()
///
double *normal_swimmer( double *n,bc WALL,smono *atom,int dimension ) {
	int i;

	if( WALL.PLANAR || ( feq(WALL.P[0],1.0) && feq(WALL.P[1],1.0) && feq(WALL.P[2],1.0) && feq(WALL.P[3],1.0) )) {
		for( i=0; i<dimension; i++ ) n[i] = WALL.A[i];
	}
	else if(( DIM == 2 && feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0)) || (DIM > 2 && feq(WALL.P[0],2.0) && feq(WALL.P[1],2.0) && feq(WALL.P[2],2.0))) {
		for( i=0; i<dimension; i++ ) n[i] = 2.*WALL.A[i]*(atom->Q[i]-WALL.Q[i]);
	}
	else for( i=0; i<dimension; i++ ) n[i] = (WALL.P[i]) *smrtPow( WALL.A[i]*(atom->Q[i]-WALL.Q[i]) , WALL.P[i]-1.0 );
	return n;
}

///
/// @brief	This subroutine applies the BC transformation to the velocity of the MD particle,
///			swimmer's version.
///
///	It transforms the velocity (the normal and tangential components) of the MD particle considering the 
/// conditions at the surface of the boundary. For instance, it can \b conserve the \b energy/momentum/
/// \b angular \b momentum using the impulse method or it can apply the rule method such as \b 
/// bounceback or \b reflection or \b periodic which does NOT necesarily conserve momentum. The BCs 
/// global variables defined in the definition.h set how particle velosity will be updated.
///
/// @param atom		The MD particle, being either the head or the middle monomer of the swimmer.
/// @param WALL 	One of the walls of the BCs that the particle is interacting with.
/// @param SS 		It specifies the type (features) of the swimmer to which this monomer belongs.
/// @param n		The normal vector to the surface of the wall.
/// @see			velBC()
///
void velBC_swimmer( smono *atom,bc *WALL,specSwimmer SS,double n[_3D] ) {
	double V[_3D],VN[_3D],VT[_3D],VR[_3D],R[_3D],zip[_3D];
	double atom_POS[_3D],atom_VEL[_3D];
	double IIpart[_3D][_3D],IIwall[_3D][_3D];
	double IMpart,IMwall,mass;	//Inverse mass
	double J=1.0;				//Impulse
	int i,j;

	//Zero everything
	for( i=0; i<_3D; i++ ) {
		V[i] = 0.0;
		VN[i] = 0.0;
		VT[i] = 0.0;
		VR[i] = 0.0;
		R[i] = 0.0;
		zip[i] = 0.0;
		atom_POS[i] = atom->Q[i];
		atom_VEL[i] = atom->V[i];
	}
	if( atom->HorM ) mass = SS.middM;
	else mass = SS.headM;

	//Calculate rotational velocity of BC
	for( i=0; i<DIM; i++ ) R[i] = atom_POS[i] - WALL->Q[i];
	crossprod( WALL->L,R,VR );
	//Begin Determining the direction of the impulse
	for( i=0; i<DIM; i++ ) V[i] = atom_VEL[i] - WALL->V[i] - VR[i];
	//Calculate normal and tangential components of velocity
	proj( V,n,VN,DIM );
	tang( V,VN,VT,DIM );

	//Point particles do not have angular momentum
	for( i=0; i<DIM; i++ ) for( j=0; j<DIM; j++ ) IIpart[i][j] = 0.;
	if( atom->HorM ) IMpart = SS.imiddM;
	else IMpart = SS.iheadM;
	if( !(WALL->DSPLC)) {
		for( i=0; i<_3D; i++ ) for( j=0; j<_3D; j++ ) IIwall[i][j] = 0.;
		IMwall = 0.0;
	}
	else if ( WALL->DSPLC ) {
		invert3x3(IIwall,WALL->I);
		IMwall = 1.0 / WALL->MASS;
	}
	else {
		printf( "Error: WALL.DSPLC must be 0 or 1.\n" );
		exit( 1 );
	}

	//The energy/momentum/angular momentum conserving impulse method
	if( WALL->COLL_TYPE == BC_IMP ) {
		//Transform the normal and tangential components
		for( i=0; i<DIM; i++ ) {
			VN[i] *= WALL->MVN;
			VN[i] += WALL->DVN*n[i];
			VT[i] *= WALL->MVT;
			VT[i] += WALL->DVT;
		}
		//Combine normal and tangential components
		V[0] = VN[0] + VT[0];
		V[1] = VN[1] + VT[1];
		if( DIM >= _3D ) V[2] = VN[2] + VT[2];
		//The impulse's direction is the difference between the two (final minus initial)
		for( i=0; i<DIM; i++ ) V[i] -= ( atom_VEL[i] - WALL->V[i] - VR[i] );
		//Normalize V. We only want the unit vector. Conservation will determine its magnitude
		norm( V,DIM );

		J = impulse( V,atom_VEL,WALL->V,atom_POS,WALL->Q,zip,WALL->L,IMpart,IMwall,IIpart,IIwall,atom_POS,WALL->E );
	}
	//The rule method such as bounceback or reflection or periodic which does NOT necesarily conserve momentum
	else if( WALL->COLL_TYPE == BC_SURF_RULES ) {
		//Transform the normal and tangential components
		for( i=0; i<DIM; i++ ) {
			VN[i] *= WALL->MVN;
			VN[i] += WALL->DVN*n[i];
			VT[i] *= WALL->MVT;
			VT[i] += WALL->DVT;
		}
		//Combine normal and tangential components
		V[0] = VN[0] + VT[0];
		V[1] = VN[1] + VT[1];
		if( DIM >= _3D ) V[2] = VN[2] + VT[2];
		//Move the velocity out of the particle's rest frame and back into the lab frame
		for( i=0; i<DIM; i++ ) V[i] += WALL->V[i] + VR[i];
		//Make V the change in momentum instead of the velocity by subtracting the old velocity and times mass
		for( i=0; i<DIM; i++ ) {
			V[i] -= atom_VEL[i];
			V[i] *= mass;
		}
		//Set J=1 because the total change in momentum (i.e. impulse) is in V.
		J = 1.0;
	}
	else if( WALL->COLL_TYPE == BC_THERMO_SURF ) {
		//Transform the normal and tangential components
		printf("Error: Thermalized wall not implemented for swimmers.\n");
		exit(1);
	}
	else if( WALL->COLL_TYPE == BC_HALF_RULES ) {
		printf("Error: Rule-based solvent/BC collisions at t/2 not yet operational.\n");
		exit(1);
	}
	else if( WALL->COLL_TYPE == BC_THERMO_HALF ){
		printf("Error: Probabilistic solvent/BC collisions at t/2 not yet operational.\n");
		exit(1);
	}
	else {
		printf("Error: Solvent/BC Collision type unrecognized.\n");
		exit(1);
	}

	//Use the impulse to set the velocity of particleMPC
	for( i=0; i<_3D; i++ ) atom->V[i] += V[i] * J * IMpart;
	if( WALL->DSPLC ) {
		//Set the velocity of BC
		for( i=0; i<DIM; i++ ) WALL->V[i] -= V[i] * J * IMwall;
		//Set the angular velocity of BC
		//Since VN isn't being used, use VN as the difference
		for( i=0; i<_3D; i++ ) VN[i] = atom_POS[i] - WALL->Q[i];
		//Since VT isn't being used, use VT as the crossprod.
		crossprod( VN,V,VT );
		//dotprodmat( VT,IIwall,VN,_3D );
		dotprodMatVec( IIwall,VT,VN,_3D );
		for( i=0; i<_3D; i++) WALL->dL[i] -= VN[i] * J;
	}
}

///
/// @brief	This subroutine applies the BC transformation to the position of the MD particle
/// 		swimmers' version.
///
/// It updates the position of the MD particle by applying the normal and tangential displacements 
/// specified in the input file, when crossing a periodic boundary.
///
/// @param atom 	The MD particle, being either the head or the middle monomer of the swimmer. 
/// @param WALL 	One of the walls of the BCs that the particle is interacting with.
/// @param n 		The normal vector to the surface of the wall.
/// @note			It does NOT stream! That is done in a seperate routine.
/// @see			posBC()
///
void posBC_swimmer( smono *atom,bc WALL,double n[_3D] ) {
	double PN[_3D],PT[_3D],temp[_3D];
	int i;

	for( i=0; i<DIM; i++ ) temp[i] = atom->Q[i];
	proj( temp,n,PN,DIM );		//Calculate normal component of the position
	tang( temp,PN,PT,DIM );		//Calculate tangential component of position
	//Transform the position
	for( i=0; i<DIM; i++ ) {
		PN[i] += WALL.DN * n[i];//Transform normal component
		PT[i] += WALL.DT;	//Transform tangential component
	}
	//Combine normal and tangential components
	for( i=0; i<DIM; i++ ) atom->Q[i] = PN[i] + PT[i];
}
