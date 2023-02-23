#ifndef MTOOLS_H
#define MTOOLS_H

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ DECLARE FUNCTIONS *********** */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/*
   These functions are nonspecific mathematical
   tools that are used through out the calculations.
*/
//Floats
int feq(double x,double y);
int fneq(double x,double y);

// Vectors
double length( double x[],int dimension );
void norm( double x[],int dimension );
void normCopy( double xin[],double xout[],int dimension );
void normalplane( double x[3], double y[3], double n[3] );
double dotprod( double x[], double y[],int dimension );
void dotprodMatVec( double M[][3],double v[],double result[],int dimension );
void dotprodVecMat( double v[], double M[][3],double result[],int dimension );
void dotprodMatMat( double A[][3],double B[][3],double result[][3],int dimension );
void crossprod( double x[3], double y[3], double result[3] );
void oldcrossprod( double x[3], double y[3], double result[3] );
void proj( double v[],double n[],double VP[],int dimension );
void tang( double v[],double VN[],double VT[],int dimension );
double cosang( double v1[],double v2[],int dimension );
double atan2( double y,double x );
double absAngle( double v1[], double v2[], int dimension );
double signedAngle( double v1[], double v2[], int dimension );

// Tensors
void outerprod( double x[], double y[], double result[][3],int dimension );
int levicivita( int i, int j, int k );
double det2x2( double m[2][2] );
double det3x3( double m[3][3] );
double determinant( double **a,int n );
double trace( double **a,int n );
void invert2x2( double m[2][2] );
double cofactor3x3( double m[3][3],int i,int j );
void invert3x3( double m_inv[3][3],double m[3][3] );
void parallelaxis( double I[][3],double R[],double M,int dimension );

// Distances
double distpoints( double P1[3],double P2[3],int dimension );
double distsurf( bc WALL,double P[3] );
double distplane( bc WALL,double x, double y, double z );
double pythag( double x, double y );
double SIGN( double x,double y );

// Frame of reference with respect to a moving BC or the system in general
void restframe( double V[],bc WALL,int dimension );
void labframe( double V[],bc WALL,int dimension );
void galileantrans( particleMPC *pp,bc WALL[],simptr simMD,spec SP[],double KBT,double VEL[],int POP,int NBC,int MDmode,int dimension );
void zeroExtraDims( particleMPC *pp,bc WALL[],simptr simMD,int GPOP,int NBC,int MDmode,int dimension );

// Properties of BCs
void mominert( bc *body,int XYZ[],int dimension );
void dim_vol( bc *body,int XYZ[],int dimension );
double latticeEstVol( bc *body,int XYZ[],int dimension );
void latticeEstMomInert( bc *body,int XYZ[],int dimension );

// Binning
void histbin( double values[],int hist[BINS],double minRange,double maxRange,int POP );
double stdNum( cell ***CL,int GPOP,int XYZ[3],int XYZ_P1[3] );

// Conserved Values
void conservation( double VA[],int MA,double QA[],double VB[],int MB,double QB[],double WB[],double IB[3][3],int dimension );

//Surfaces
double surf_func( bc WALL, double POS[],int dimension );
double non4foldSymmCalcW( bc WALL,double POS[], int dimension );

// Eigenvalues and vectors
void eigenvalues2x2( double **m,double eigval[] );
int eigenvectors2x2( double **m,double eigval[],double eigvec[][2] );
void eigenvalues3x3( double **m,double eigval[] );
void eigenvectors3x3( double **m,double eigval[],double eigvec[][3] );
int solveEigensystem( double **m,int dimension,double eigval[] );

//Derivatives
double centredDeriv( double xM1,double xP1,double dt );
double forwardDeriv( double x0,double xP1,double dt );
double backwardDeriv( double x0,double xM1,double dt );
double simps( double F[],double dx,int n );

//Rotations
void skewSymmetricCrossProductMatrix( double *v,double result[][3] );
void rotationMatrix( double rotMat[][3],double vx[][3],double c,double s );
void findRotationMatrix( double rotMat[][3],double *original,double *final );
void rodriguesRotation( double vec[],double rotAx[],double theta );
void setRotMatrix3D( double M[][3],double angx,double angy,double angz );
void setRotMatrix2D( double M[][3],double angz );

//Spatial correlation functions
void dirdirCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension );
void densdensCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension );
void orderorderCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension );
// void phiphiCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension );
void velvelNormedCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension );
void vortvortNormedCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension );
void velvelCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension );
void vortvortCorr( cell ***CL,int maxXYZ,int XYZ[3],double *avCorr,int dimension );
void normCorr( double *corr,int maxXYZ );
void FTspectrum( double *corr,double *spect,int maxXYZ,int dimension );

//Checks
int checkNAN_vec( double vec[],int dimension );
void checkNAN_Q( cell ***CL,int XYZ_P1[3],int pauseFlag,int dimension );
void checkNAN_V( cell ***CL,int XYZ_P1[3],int pauseFlag,int dimension );

//Optimisation related functions
double smrtPow(double x, double y);
#endif
