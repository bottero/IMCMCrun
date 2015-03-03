/*
 * structures.h
 *
 *  Created on: 27 nov. 2014
 *      Author: abottero (alexis DOT bottero aT gmail DOT com)
 */

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "defines.h"
#include "tab3d.hpp"

typedef struct{       // To store a Markov chain state
  std::vector<double> params;                // To store the parameters that describe the state
  double E;                                  // To store the energy E of the state
}State;

typedef struct{       // To store a Markov chain
  std::vector<State> states;
  int i;                                    // To store the number of the chain
  double T;                                 // To store the temperature of the chain
  int nc;                                   // Number of components that will be modified at each iteration
  int at;                                   // Number of accepted MH transitions
  int rt;                                   // Number of rejected MH transitions
  int od;                                   // Number of "out of domain"s
  int ps;                                   // Number of proposed swaps
  int as;                                   // Number of accepted swaps
  int rs;                                   // Number of rejected swaps
  std::vector<double> deltaParameters;      // Maximum amplitude of variation for each parameter
  std::vector<double> accProba;             // Instantaneous acceptance probability
  std::vector<double> swapProba;            // Instantaneous swapping probability
  std::vector<std::vector<double> > profilesP;   // To keep the profiles and compute the quantiles 
  std::vector<std::vector<double> > profilesS;   // To keep the profiles and compute the quantiles
  std::vector<double> maxP;                 // To keep maximum P wave velocity values investigated by the chain
  std::vector<double> maxS;                 // To keep maximum S wave velocity values investigated by the chain
  std::vector<double> minP;                 // To keep minimum P wave velocity values investigated by the chain
  std::vector<double> minS;                 // To keep minimum S wave velocity values investigated by the chain
  // Ex : profilesP[15] is a vector of size "n" that contains the history of the P wave velocity at the point 15 of the profile
  std::vector<double> averageP, varP, qInfP, qSupP; // To save the medium value of P wave velocity, the variance and the quantiles
  std::vector<double> averageS, varS, qInfS, qSupS; // To save the medium value of S wave velocity, the variance and the quantiles
}Chain;

typedef struct{       // To store the features of a swap
  int idxChain;                             // Index of the chain which has swapped (0 : lowest temperature chain...)
  int idxState;                             // Index of the state of the chain who has been swapped
  int sp;                                   // Index of the state of the higher temperature chain who has been chosen for the swap
}SwapFeatures;

typedef struct{       // To store parallel Markov chains
  std::vector<Chain*> chains;                 // To store the Markov chains
  std::vector<std::vector<long double> > SCI; // To store the matrix of importance weights
  std::vector<SwapFeatures> swapHist;         // To store the history of the swaps
  std::vector<double> maxP;                 // To keep maximum P wave velocity values investigated by the algorithm
  std::vector<double> maxS;                 // To keep maximum S wave velocity values investigated by the algorithm
  std::vector<double> minP;                 // To keep minimum P wave velocity values investigated by the algorithm
  std::vector<double> minS;                 // To keep minimum S wave velocity values investigated by the algorithm
  std::vector<int> idxE;  // To keep indices of best models
  std::vector<double> bestE;  // To keep energies of best models
  std::vector<int> chainBestE;  // To keep best models chains
}Run;

typedef struct{       // To store coordinates expressed in meters
  int index;
  int chainNumber;
  double E;
}IndexBestEnergies;

typedef struct{       // To store coordinates expressed in meters
  double x;
  double y;
  double z;
}Coordinate;

typedef struct{       // To store the arrival times at a receiver
  std::vector<double> timesP;               // First arrival times P at each receiver : 1-32 shot 1, 33-65 shot 2 ...
  std::vector<double> timesS;               // First arrival times S at each receiver : 1-32 shot 1, 33-65 shot 2 ...
}ArrivalTimes;

typedef struct{       // 3D arrays to store the velocity model (P and S waves velocities at each point)
  tab3d<double>* velP;
  tab3d<double>* velS;
  int nx;                                   // (Number of cell in direction x) + 1 -> number of borders
  int ny;                                   // (Number of cell in direction y) + 1
  int nz;                                   // (Number of cell in direction z) + 1
  float dx;                                 // Interval in direction x
  float dy;                                 // Interval in direction y
  float dz;                                 // Interval in direction z
  double xmin;                              // minimum of x coordinate
  double ymin;                              // minimum of y coordinate
  double zmin;                              // minimum of z coordinate
}VelocityModel;

typedef struct{       // To store all that is related with MPI
  int rc;
  int len;
  int nb_process;                        // Number of processes
  int rank;                              // Rank of the process considered
  #ifdef PAR  
  char hostname[MPI_MAX_PROCESSOR_NAME];
  MPI_Status status;
  #endif
}MpiConfig;

typedef struct{       // To store data characteristics
  ArrivalTimes times;                    // Arrival times at each station for each shot : 1-32 shot 1, 33-65 shot 2 ...
  std::vector<Coordinate> coordStations; // Coordinates of the stations
  std::vector<Coordinate> coordShots;    // Coordinates of the shots
  std::vector<double> firstGuessP;       // FirstGuess for P waves (size nz-1)
  std::vector<double> firstGuessS;       // FirstGuess for S waves (size nz-1)
  std::vector<double> filtFirstGuessP;          // Filtered First guess for P waves velocity (size nz-1)
  std::vector<double> filtFirstGuessS;          // Filtered First guess for S waves velocity (size nz-1)
  std::vector<double> zpFirstGuess;             // Corresponding depths (size nz-1)
  std::vector<double> realP;             // If ANALYTICAL_RUN == 1 store the real P profile (size nz-1)
  std::vector<double> realS;             // Idem for real S profile (size nz-1)
  std::vector<double> zp;                // Corresponding depths (size nz-1) (we verify that zp=zpFirstGuess)
  std::vector<double> z;                 // Corresponding border depths (size nz)
  std::vector<double> zFilt;             // Depths for filtered profiles (we need less points to describe the filtered profiles properly) (size nzFilt)
  std::vector<double> zFiltp;            // Corresponding depths (size nzFilt-1)
  std::vector<double> minParameters;     // Minimum for each parameter
  std::vector<double> maxParameters;     // Maximum for each parameter
  std::vector<int> indexParameters;      // The params[i] refer to the wavelet coefficient number indexParameters[i] of the decomposition of the profiles 
  std::vector<int> staticParameters;     // Indexes of the non varying parameters : Some parameters can be fixed to a value to reduce the dimension
  double Ep;                             // Energy of the prior
  double sigmaP;                         // Standard deviation of the P arrival times picked
  double sigmaS;                         // Standard deviation of the S arrival times picked
  int nx;                                // (Number of cell in direction x) + 1
  int ny;                                // (Number of cell in direction y) + 1
  int nz;                                // (Number of cell in direction z) + 1 for the first guesses (and the real profiles)
  int nzFilt;                            // (Number of cell in direction z) + 1 for the filtered profiles < nz (we need less 
                                         // points to describe the filtered profiles properly)
  float dx;                              // Interval in direction x
  float dy;                              // Interval in direction y
  float dz;                              // Interval in direction z for the first guesses (and the real profiles)
  float dzFilt;                          // Interval in direction z for the filtered profiles > dz (we need less 
                                         // points to describe the filtered profiles properly)
  double xminGrid;                       // minimum of x coordinate
  double yminGrid;                       // minimum of y coordinate
  double zminGrid;                       // minimum of z coordinate
}Data;

typedef struct{       // To store the configuration of the run
  time_t startTime;           // To store the starting time of the program
  MpiConfig mpiConfig;        // To store MPI's stuffs
  Data data;                  // Data used
  int npu;                    // Number of wavelet coeffs used to describe P wave velocity (and S wave velocity)
  int nbt ;                   // Number of temperatures (Number of Markov chains in parallel)
  int nit;                    // Number of iterations
  int nSweeps;                // Number of sweeps in the Eikonal
  float epsin;                // For the Eikonal, radius in number of grid points arround source where spherical approximation will be used
  double di,df;               // Control the amplitude of the steps (At T=Tmax deltaState=L/DI and at T=1 deltaState=L/DF)
  double tmax ;               // Temperature max
  double pee;                 // Probability of allowing swapping
  int seed;                   // Seed of the random numbers generators
  std::vector<int> nc;
  // Number of component to modify in a MH iteration as a function of T. Example with 100 parameters and 4 temperatures :
  // nc={5,20,50,80}; -> at each iteration 5 parameters of the lowest temperature chain can change (80 parameters of the highest).
  std::vector<double> T;                    // Temperature ladder (T[0]=1)
  double A;                              // Factor used to create the prior from the wavelet transform of the first guess
                                         // If coeffs[i] > 0 :
                                         // (1-A)coeffs[i] < wavelet coefficients number i generated by the algorithm < (1+A)coeffs[i]
  int nPriorProfiles;                    // Number of profiles from a priori space generated in initialization
  std::string code;                         // To store a combination used to distinguish files from different runs
  std::string filesDir;                     // Directory where the data files are stored
  std::string outputDir;                    // Path to the output directory
  std::string confFile;                     // Path to the configuration file (ex : /home/abottero/charon.conf)
  std::string name_of_real_profile_P;
  std::string name_of_real_profile_S;
  std::string name_of_times_file;
  std::string name_of_stations_file;
  std::string name_of_shots_file;
  std::string name_of_prior_features_file;
  std::string name_of_first_guess_P_file;
  std::string name_of_first_guess_S_file;
  /************Quantile Options***************/
  double qp;                                // Ratio of the values on the quantile chosen (ex : 0.95 -> 95%)
  /************Wavelet variables**************/
  std::string wavelet;                           // Wavelet used for the wavelet based filering
  int ndwts;                                    // Number of DWT stages for the wavelet transform
  std::vector<double> coeffsP, flagP;       // For the wavelet transform (coeffsP is used to store temporary coefficients)
  std::vector<int> lengthP;                 // For the wavelet transform
  std::vector<double> coeffsS, flagS;       // For the wavelet transform
  std::vector<int> lengthS;                 // For the wavelet transform
  /************Other options***************/
  int swaves;
  int verbose1;
  int verbose2;
  double coordTol;
  int buildPrior;
  int findOptimumGrid;
  int shotNumberRef;
  int nxref;
  int nyref;
  int numberOfXpointsToDescribeTheSmallestDistanceDefault;
  int numberOfYpointsToDescribeTheSmallestDistanceDefault;
  int nzfiltDefault;
  int analyticalRun;
  int recalculateT0;
  int useDefaultSeed;
  int defaultSeed;
  int test;
  double minEtest;
  double maxEtest;
  int view;
  std::vector<int> noXptdtsdVec;
  std::vector<int> noYptdtsdVec;
  std::vector<int> nzFiltVec;
  int nBestProfiles;
  int computeResiduals;
  int iterationsResiduals;
  int iterationsBestProfiles;

}Configuration;

/*
struct dbl_cmp {
    dbl_cmp(double v, double d) val(v), delta(d) { }
    inline bool operator()(const double &x) const {
        return fabs(x-val) < delta;
    }
private:
    double val, delta;
};*/

#endif /* STRUCTURES_H_ */
