/*
 * defines.h
 *
 *  Created on: 27 nov. 2014
 *      Author: abottero (alexis DOT bottero aT gmail DOT com)
 */
 
#ifndef DEFINES_H_
#define DEFINES_H_

//********** You should not have to change anything on that file **********

#define TINYVAL 1E-8          // To test if 2 values are different we do : fabs(a-b) > TINYVAL. Coordinates must be >> TINYVAL
#define TMIN 1                // Min temperature
#define NAME_OF_CONFIGURATION_FILE "config.cfg"
#define TRESH 50              // Max and min values are computed from iteration TRESH

extern "C" {void fteik_(double *vel, double *times,int *nz,int *nx,int *ny,float *zsin,float *xsin,float *ysin,float *dzin,float *dxin,float *dyin,int *nsweep,float *epsin);}
                              // Fortran subroutine for the Eikonal
                                    
#ifdef PAR  // PAR is defined during the Makefile call : make PAR=yes. See the Makefile
#include "mpi.h"
#endif

#endif /* DEFINES_H_ */
