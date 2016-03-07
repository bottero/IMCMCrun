//===================================================================================================
// Name        : RealisticDataMCMC.cpp
// Author      : Alexis Bottero (alexis DOT bottero aT gmail DOT com)
// Description : Interacting Monte-Carlo Markov chains applied to a realistic dataset
//===================================================================================================

#include <iostream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "rngs.h"     // Libraries with random functions
#include "rvgs.h"     //               "
#include "functions.h"
#include "defines.h"
#include "structures.h"
#include "initializations.h"
#include "generalFunctions.h"
#include "filesAndControl.h"

int main( int argc, char *argv[]) // argc : Number of arguments, argv : the arguments. !! argv[0] is the name of the program !!
{
  //*************************************************  Initialization : *********************************************************//
  
  Configuration config;                     // Initialize the configuration structure that will contain the config chosen for the run
  config = parameters_setting(argc,argv);   // Store the configuration of the run in the structure "config"
  loadData(&config);                        // Load the data (arrival times, sources and receivers positions) contained on the files
  designSmallGrid(&config);                 // Calculate optimum nx,ny,nzFilt,dx,dy,dzFilt to reduce computation time
  if (config.calculateTimesForFirstGuess) { // With this option we just run the eikonal for the profile(s) given as first guess(es)
    calculateTimesForFirstGuess(&config);   // ...
    MPI_Finalize();                         // ...
    exit(0);                                // ... and we terminate!
  }
  if (config.buildPrior)                  //
    buildPrior(&config);                  // Build the prior features
  generate_profiles_from_prior(&config);  // Generate config.nPriorProfiles profiles from a priori space in config.outputDir/priorCurvesXXX/
  if (config.analyticalRun)               // If we perform an analytical run :
    createDataset(&config);               // ... create the arrival-times from the real profile
  Run run;                                // Create the run, it will contain all the chains
  run = init_run(&config);                // Initialize it (fill it with chains of initialized states)
  writeStatus(&run,&config);              // Create the data files (chainI.XXX.dat) and write the first line on it.
  write_config(&run,&config);             // Create a file config.XXX.dat and write the configuration parameters on it

  //**********************************************  End of initialization : *****************************************************//

  for (int n=0;n < config.nit;n++) { // Loop on the iterations
    printEvolution(n,&run,&config);
    for (int i=config.nbt-1;i>=0;i--) { // Loop on the different temperatures chains
      if (i == config.nbt-1) {              // If we consider the highest T chain...
       // priorIteration(run.chains[i],&config);              // (...we perform an iteration on the prior.)
        iterationMHindependent(run.chains[i],&config);      // ...we perform an independent sampler transition.
      }
      else {                                // Otherwise, if we consider the low temperature chains...
        double p=Random();
        if (p < (1-config.pee))                             // _with probability 1-pee we perform a classical MH transition.
          iterationMH(run.chains[i],&config);               // Iteration Metropolis Hasting
        else                                                // _with probability pee we suggest an IR-swap.
          importanceSamplingSwap(&run,i,&config);           // ImportanceSamplingSwap
      }
      // updateAverageProfiles(run.chains[i],&config,n+1);  // Old version: Update the average, variance and quantiles profiles (n+1 because we have to count the profile created during initialization)
    }
    updateAverageProfiles(&run,&config,n+1);  // Update the average, variance and quantiles profiles (n+1 because we have to count the profile created during initialization)
    updateSCI(&run,n);          // Add a new line to SCI (importance weights (cumulative sum) and normalization coefficients)
    writeFiles(&run,&config,n); // Write informations on files
  }
  summary(&run,&config);       // Display a summary of the run
  write_summary(&run,&config); // Write a summary of the run at the end of the file config.XXX.dat
  printTime(&config);          // Print the final time
  finalizeRun(&run);           // Free memory allocated for the chains and finalize MPI

  //************************************************  End of the program ********************************************************//

  return 0;
}
  
