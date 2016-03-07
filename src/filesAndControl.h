#ifndef FILESANDCONTROL_H_
#define FILESANDCONTROL_H_

#include "structures.h"

void checkConfig(const Configuration* config);
// Check if the configuration is correct. TODO make this function
void write_config(const Run* run, const Configuration* config);
// Create a file config.XXX.dat and write the configuration parameters in it.
void copyDataFiles(const Configuration* config);
// Copy the data files used during the simulation in OUTPUT_FILES
void printEvolution(const int n, const Run* run,const Configuration* config);
// Print the iteration number n. TODO : add a % done and at when the run will finish.
void writeFiles(Run* run, Configuration* config, int n);
// Write informations on files
void writeStatus(Run* run,Configuration* config);
// If it does not exist, creates the data files (chainI.XXX.dat). Writes a new line on it.
void writeSwap(Run* run,Configuration* config);
// If it does not exist, creates the data file (ll.XXX.dat). Writes a new line on it. It store the swapping history
void writeSCI(Run* run,Configuration* config);
// Print the Importance Weights matrix SCI in sci.XXX.dat
void writeStats(Run* run,Configuration* config);
// Write the statistics of the chains (probability of swapping...) in statsI.XXX.dat
// In the order : at rt od ps as rs acceptance probability swapping probability : see structure Chain or README
void printStatus(Run* run,const Configuration* config,int view);
// Print the status of the run, display the "view" last states
void summary(Run* run,const Configuration* config);
// Display a summary of the run
void write_summary(const Run* run,const Configuration* config);
// Write a summary of the run at the end of the file config.XXX.dat
void writeAverages(const Run* run, const Configuration* config);       
// Write the average, variance and quantiles profiles on files
void writeBestProfiles(Run* run, const Configuration* config);
// Write the best profiles on files
void computeResidualsForBestModel(Run* run, Configuration* config);
// Recompute the eikonal (in serial) for the best model and save the residuals on a file
void writeMinMaxProfiles(Run* run, const Configuration* config); 
// Update min and max velocities explored by the run, write them on file
void write_one_column_file(const std::vector<double>* column, const std::string name_of_file);
// Write a a vector into a column file
void write_two_columns_file(const std::vector<double>* column1, const std::vector<double>* column2, const std::string name_of_file);
// Write a two columns file from two vectors of the same size
void write_three_columns_file(const std::vector<int>* column1, const std::vector<double>* column2, const std::vector<double>* column3, const std::string name_of_file);
// Write a three columns file from three vectors of the same size -> just used for prior features
void write_four_columns_file(const std::vector<int>* column1, const std::vector<double>* column2, const std::vector<double>* column3, const std::vector<double>* column4, const std::string name_of_file);
// Write a four columns file from four vectors of the same size

#endif /* FILESANDCONTROL_H_ */
