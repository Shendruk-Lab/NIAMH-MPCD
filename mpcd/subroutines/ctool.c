///
/// @file
/// @brief Helpful routines
///

# include<stdio.h>
# include<math.h>
# include <stdarg.h>

/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
/* ************ HELPFUL ROUTINES ************ */
/* ****************************************** */
/* ****************************************** */
/* ****************************************** */
///
/// @brief		Waits for the user to press enter.
///		
void wait4u() {
	printf( "\nPress enter to continue\n" );
	while ( 1 ) {
		if( '\n' == getchar() )
		break;
	}
}

///
/// @brief			    Prints a warning message if the value is within the precision cutoff for zero. 
///
/// @param test			The value that will be checked.
/// @param zero 		The check value for zero.
/// @param mult 		The coefficient that changes the magnitude of `zero`.
///
void zerowarning( double test, double zero,double mult ) {
	if( fabs( test ) <= zero*mult ) printf( "Warning: Value %e is within %lf times magnitude of zero check value of %e.\n",test,mult,zero );
}

// NaN "hunting" routines. Primarily for debugging (set breaks here!).

///
/// @brief Uses the isnan() preprocessor macro to check for NaNs.
///
/// A wrapper for the isnan() preprocessor macro. This is particularly useful to use as a break point in a debugger.
///
/// @param x Check this double for NaNs.
/// @return 1 if NaN, 0 if not.
///
int isNaN(double x) {
    if (isnan(x)) {
        printf("NaN detected.\n");
        return 1;
    }
    else return 0;
}

///
/// @brief Uses the isnan() preprocessor macro to check for NaNs on an array.
///
/// A wrapper for the isnan() preprocessor macro. This is particularly useful to use as a break point in a debugger.
///
/// @param x The double array to check for NaNs.
/// @param n The length of the array.
/// @return 1 if NaN, 0 if not.
/// @see isNaN(double x)
///
int isNaNs(double *x, int n) {
    for (int i=0; i<n; i++) {
        if (isnan(x[i])) {
            printf("NaN detected at index %d.\n",i);
            return 1;
        }
    }
    return 0;
}

///
/// @brief Function that zeros any vector.
///
/// This function sets the component of the receiving vector to 0.
///
/// @param VEC The vector whose components will be zeroed.
/// @param dimension The dimension of VEC.
///
void zerovec( double VEC[],int dimension ) {
    int i;
    for( i=0; i<dimension; i++ ) VEC[i]=0.0;
}

///
/// @brief Variadic version of zerovec
///
/// This function will set all vectors following the `dim` param to zero.
/// Ex: `zerovec_v(3, _3D, vec1, vec2, vec3);` will set the 3 3D vectors to zero.
///
/// @param count The number of vectors you are zero'ing.
/// @param dim The dimension of the vectors.
/// @param ... All vectors to be zeroed, cast as `double *`, and seperated by commas.
/// @see zerovec
///
void zerovec_v(int count, int dim, ...) {
    int i;

    va_list ap;
    va_start(ap, dim);
    for (i=0; i<count; i++) {
        double *vec = va_arg(ap, double *);
        zerovec(vec, dim);
    }
    va_end(ap);
}
