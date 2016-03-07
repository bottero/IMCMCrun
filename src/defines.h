/*
 * defines.h
 *
 *  Created on: 27 nov. 2014
 *      Author: abottero (alexis DOT bottero aT gmail DOT com)
 */
 
#ifndef DEFINES_H_
#define DEFINES_H_

//********** You should not have to change anything on that file **********

#define TINYVAL 1E-8           // To test if 2 values are different we do : fabs(a-b) > TINYVAL. Coordinates must be >> TINYVAL
#define TMIN 1                 // Min temperature
#define NAME_OF_CONFIGURATION_FILE "config.cfg"
#define TRESH 50               // Max and min values are computed from iteration TRESH
#define MAX_DIFF_TO_BE_OK 2 // Tolerance to determine if we are precise enough with the eikonal. We set the velocity to 1 and we check that
// we obtain the correct distance between sources and receivers.
#define PREC 128               // Precision for exponential comparisons
#define RND MPFR_RNDU          // Rounding
/*U round toward plus infinity
D round toward minus infinity
Y round away from zero
Z round toward zero
N round to nearest (with ties to even)
*/
//#define NAN 1.0/0.0;

extern "C" {void fteik_(double *vel, double *times,int *nz,int *nx,int *ny,float *zsin,float *xsin,float *ysin,float *dzin,float *dxin,float *dyin,int *nsweep,float *epsin);}
                              // Fortran subroutine for the Eikonal

// Wavelet used for the inversion:
#define NUMBER_OF_WAVELETS 5
#define LIST_OF_WAVELETS "db8","db12","haar","coif3","coif5"
//"db6","db8","db12","haar","bior3.5","bior3.7","bior3.9","coif2","coif3","coif5"

//"haar","db1","db2","db3","db4","db5","db6","db7","db8","db9","db10","db11","db12","db13","db14", "db15","bior1.1","bio1.3","bior1.5","bior2.2","bior2.4","bior2.6","bior2.8","bior3.1","bior3.3","bior3.5","bior3.7","bior3.9","bior4.4","bior5.5","bior6.8","coif1","coif2","coif3","coif4","coif5"

#ifdef PAR  // PAR is defined during the Makefile call : make PAR=yes. See the Makefile
#include "mpi.h"
#endif

#endif /* DEFINES_H_ */
