/*
 * initialisations.h
 *
 *  Created on: 21 jan 2014
 *      Author: abottero
 */

#ifndef INITIALISATIONS_H_
#define INITIALISATIONS_H_

#include "structures.h"

Configuration parameters_setting(int argc,char *argv[]);
// Initialize the configuration structure containing all the config of the run
void readConfFile(std::string nameOfConfFile, Configuration* config);
// Read the informations contained in the configuration file
void build_temperature_ladders(Configuration* config);
// Build the temperature ladders in the configuration structure (T[0]=1)
void loadData(Configuration* config);
// Load the data (arrival times, sources and receivers positions) contained on the files
void loadArrivalTimes(Configuration* config);
// Load the first arrival times for P and S waves contained on the file tObs4096.txt
void loadStationsPositions(Configuration* config);
// Load the receivers positions contained on the file coordStats.txt
void loadShotsPositions(Configuration* config);
// Load the sources positions contained on the file coordShots.txt
void loadPriorFeatures(Configuration* config);
// Load parameters features (min, max, amplitude of variation) contained on the file priorFeatures.txt
void findOptimumFirstGuessParameterization(Configuration* config);
// From the first guess velocity profile(s) and given the parameterization scheme (store in config), determine the best set of
// coefficients describing it. Store it in config->coeffsP (and in config->coeffsS). Then store the corresponding parameterized profile
// in config->data.filtFirstGuessP (and in config->data.filtFirstGuessS).
void loadFirstGuesses(Configuration* config);
// Load first guesses on files falseFirstGuessP4096.txt, falseFirstGuessS4096.txt. Perform the discrete wavelet transform and store the results in config
// Must be done after loading priorFeatures.txt 
void buildPrior(Configuration* config);
// Filter out the first guess and record the most significant coefficient (config->data.indexParameters)
// define their maximum variation ranges (config->data.minParameters, config->data.maxParameters) during the rest of the algorithm
void generate_profiles_from_prior(Configuration* config);
// Generate config.nPriorProfiles profiles from a priori space in config.outputDir/priorProfilesXXX/
void designSmallGrid(Configuration* config);
// Calculate optimum nx,ny,nzFilt,dx,dy,dzFilt to reduce computation time
void calculateTimesForFirstGuess(Configuration* config);
// Runs the eikonal for first guess(es) curves
void findOptimumGrid(double xmaxGrid, double ymaxGrid, double dxmin, double dymin, Configuration* config);
// Filter a first guess curve iteratively and run the eikonal each time to determine nx,ny and nzFilt to use
// Then calculate dx,dy and dzFilt
void createDataset(Configuration* config);
// Create the arrival times from the real profile
State init_state();
// Initialize a markov chain state
Run init_run(Configuration* config);
// Initialize the run, create each chain at each temperature and the first state of each one.

#endif /* INITIALISATIONS_H_ */
