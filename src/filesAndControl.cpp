#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>      // std::setprecision
#include <string>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <errno.h>  // To detect mathematical errors
#include <algorithm>    // std::min_element, std::max_element
#include <iterator>
#include "functions.h"
#include "structures.h"
#include "defines.h"
#include "rngs.h"
#include "rvgs.h"
#include "generalFunctions.h"
#include "filesAndControl.h"

void checkConfig(const Configuration* config)
// Check if the configuration is correct. TODO continue this function
{
  if (config->verbose2) {
    for (int i=0; i<config->nbt;i++) 
      std::cout << "config->nc[" << i << "] = " << config->nc[i] << std::endl;
    std::cout << std::endl;
    for (unsigned int i=0; i<config->noXptdtsdVec.size();i++)
      std::cout << "config->noXptdtsdVec[" << i << "] = " << config->noXptdtsdVec[i] << std::endl;
    std::cout << std::endl;
    for (unsigned int i=0; i<config->noYptdtsdVec.size();i++)
      std::cout << "config->noYptdtsdVec[" << i << "] = " << config->noYptdtsdVec[i] << std::endl;
    std::cout << std::endl;
    for (unsigned int i=0; i<config->nzFiltVec.size();i++)
      std::cout << "config->nzFiltVec[" << i << "] = " << config->nzFiltVec[i] << std::endl;
    std::cout << std::endl;
    std::cout << "VIEW :" << config->view << std::endl;
    std::cout << "COORD_TOL :" << config->coordTol << std::endl;
    std::cout << "DEFAULT_SEED :" << config->defaultSeed << std::endl;
  }
}
void write_config(const Run* run, const Configuration* config)
// Create a file config.XXX.dat and write the configuration parameters in it.
{
  if (config->mpiConfig.rank == 0) {       // (In case of parallel implementation just one process has to create files)
    std::string conf(config->outputDir+"config."+config->code+".dat");  // We will obtain /home/abottero/config.321.dat for example
    std::ofstream file(conf.c_str()); // Create the file
    if(file == NULL) {
      std::cout << "IMPOSSIBLE TO CREATE config.XXX.dat" << std::endl;
      exit(0);
    }
    file << "*********************** Importance resampling MCMC ***********************" << std::endl;
    file << "*************************** REALISTIC DATASET ****************************" << std::endl;
    file << "***************** configuration file for run number "<< config->code <<" ******************" << std::endl;
    char buffer[80];
    convertTime(buffer, config->startTime);
    file << "Starting time : " << buffer << std::endl;
    file << "Seed of the random number generator = " << config->seed;
    if(config->useDefaultSeed)
      file << " (config->defaultSeed)";
    file << std::endl;
    if(config->test) {
      file << "THIS IS A TEST RUN ( we run the program without enter the Eikonal but generating random energies uniformly ";
      file << "distributed between : "<< config->minEtest << " and " << config->maxEtest << ")" << std::endl;
    }
    #ifdef PAR  
    file << "                           === MPI MODE ===" << std::endl;
    file << config->mpiConfig.nb_process << " processes running on : " << config->mpiConfig.hostname << std::endl;
    #endif
    file << std::endl;

    file << "D A T A S E T  U S E D" << std::endl << std::endl;
    if(config->swaves)
      file << "  S waves are calculated" << std::endl;
    if(config->analyticalRun) {
      file << "It is an analytical run :" << std::endl; 
      file << "For this simulation we already know the real velocity model, it has been loaded from the ";
      if(config->swaves)
        file << "files : "+config->name_of_real_profile_P+" and "+config->name_of_real_profile_S;
      else
        file << "file : "+config->name_of_real_profile_P;
      file << std::endl;
      file << "These are the 5 first lines : (z vp vs)" << std::endl;
      if(config->swaves) {
        for(int i=0;i<5;i++)
          file << config->data.z[i] << " " << config->data.realP[i] << " " << config->data.realS[i] << " " << std::endl; 
      }
      else {
        for(int i=0;i<5;i++)
          file << config->data.z[i] << " " << config->data.realP[i] << " " << std::endl; 
      }
      file << std::endl;
      file << "(dz = " << fabs(config->data.z[1]-config->data.z[0]) << ")" << std::endl << std::endl; // = (zmax-zmin)/(length(z)-1)
      file << "  The sources and receivers positions are on the files : ";
      file << config->filesDir+config->name_of_shots_file << config->filesDir+config->name_of_stations_file;
      file << std::endl;
    }
    else {
      file << "  The data used and the sources and receivers positions are on the files : ";
      file << config->filesDir+config->name_of_times_file<< ", " << config->filesDir+config->name_of_shots_file << config->filesDir+config->name_of_stations_file;
      file << std::endl;
    }
    file << "  Number of shots : " << config->data.coordShots.size() << std::endl;
    for (int i=0;i<(int)(config->data.coordShots.size());i++) {
      file << "    Shot number "<< i;
      file << " at (x=" << config->data.coordShots[i].x << ", y=" << config->data.coordShots[i].x << ", z=" << config->data.coordShots[i].z;
      file << ")" << std::endl;
    }
    file << "  Number of stations : " << config->data.coordStations.size() << std::endl;
    for (int i=0;i<(int)(config->data.coordStations.size());i++) {
      file << "    Station number "<< i;
      file << " at (x=" << config->data.coordStations[i].x << ", y=" << config->data.coordStations[i].x << ", z=" << config->data.coordStations[i].z;
      file << ")" << std::endl;
    }
    if(!config->analyticalRun) {
      file << "  Number of P waves arrival times in"+config->name_of_times_file+" : " << config->data.times.timesP.size() << std::endl;
      if(config->swaves)
        file << "  Number of S waves arrival times in the same file : " << config->data.times.timesS.size() << std::endl;
    }
    file << "  Standard deviation of P wave arrival times : " << config->data.sigmaP << std::endl;
    if(config->swaves)
      file << "  Standard deviation of S wave arrival times : " << config->data.sigmaS << std::endl;
    file << "  Number of points in P wave velocity first guess : " << config->data.firstGuessP.size() << std::endl;
    if(config->swaves)
      file << "  Number of points in S wave velocity first guessS : " << config->data.firstGuessS.size() << std::endl;
    
    file << std::endl << "P R I O R  K N O W L E D G E  O N  T H E  P A R A M E T E R S" << std::endl << std::endl;
    if (config->waveletParameterization) {
      file << std::endl << "Wavelet based parameterization" << std::endl;
      if (config->useAllWavelets) {
        file << "  Using wavelets (" << config->nWavelets << " in the list):";
        file << " ";
        for (int wavelet = 0; wavelet < config->nWavelets; wavelet++) // Loop on the wavelets used
          file << config->listOfWavelets[wavelet] << " ";
        file << std::endl;
      }
    }
    else
      file << std::endl << "Layers based parameterization" << std::endl;
    file << "  Number of parameters : " << config->data.minParameters[0].size();
    if (config->useAllWavelets and config->waveletParameterization)
      file << " plus the wavelet type.";
    std::cout << std::endl;
    file << "  (DI : " << config->di << ", DF : " << config->df << ")" << std::endl;
    for (int wavelet = 0; wavelet < config->nWavelets; wavelet++) { // Loop on the wavelets used
      if (config->useAllWavelets and config->waveletParameterization)
        file << "  Wavelet:" << config->listOfWavelets[wavelet] << std::endl;
      file << "                  min     max     ";
      for (int i=0;i<config->nbt;i++) 
        file << "delta for T" << i << "   ";
      file << std::endl;
      for (int i=0;i<(int)(config->data.minParameters[wavelet].size());i++) {
        file << "  Parameter "<< i << " :";
        file << " " << config->data.minParameters[wavelet][i] << "  " << config->data.maxParameters[wavelet][i] << "     ";
        for (int j=0;j<config->nbt;j++) 
          file << run->chains[j]->deltaParameters[wavelet][i] << "         ";
        if (config->waveletParameterization) {
          if(config->swaves && (i >= (int)((double)config->data.minParameters[wavelet].size()/2.0) )) // First params are for P waves
            file << "  (refer to wavelet coefficient number " << config->data.indexParameters[wavelet][i] << " of the S wave velocity first guess)" << std::endl;
          else
            file << "  (refer to wavelet coefficient number " << config->data.indexParameters[wavelet][i] << " of the P wave velocity first guess)" << std::endl;
        }
      }
      std::cout << std::endl;
    }
    if((int)(config->data.staticParameters.size()) != 0) {
      file << "  Static parameters :";
      for (int i=0;i<(int)(config->data.staticParameters.size());i++) 
        file << config->data.staticParameters[i] << " ";
    }
    file << std::endl;
    file << "  Energy of the prior : " << config->data.Ep << std::endl;    

    file << std::endl << "P A R A M E T E R S  O F  T H E  E I K O N A L" << std::endl << std::endl;
    file << "  nx: " << config->data.nx << "  ny: " << config->data.ny << "  nz: " << config->data.nz << std::endl;
    file << "  nzFilt: " << config->data.nzFilt << std::endl;
    file << "  dx: " << config->data.dx << "  dy: " << config->data.dy << "  dz: " << config->data.dz << std::endl;
    file << "  dzFilt: " << config->data.dzFilt << std::endl;
    file << "  xmin: " << config->data.xminGrid << "  xmax: " << config->data.xminGrid+(config->data.nx-1)*config->data.dx << std::endl;
    file << "  ymin: " << config->data.yminGrid << "  ymax: " << config->data.yminGrid+(config->data.ny-1)*config->data.dy << std::endl;
    file << "  zmin: " << config->data.zminGrid << "  zmax: " << config->data.zminGrid+(config->data.nzFilt-1)*config->data.dzFilt << std::endl;
    file << "  Number of sweeps : " << config->nSweeps << std::endl;
    file << "  EPSIN : " << config-> epsin << std::endl;

    file << std::endl << "O U T P U T  P A R A M E T E R S" << std::endl << std::endl;
    file << "  Output directory : " << config->outputDir << std::endl;
    file << "  Quantile chosen : " << config->qp << std::endl;
    
    file << std::endl << "C O N F I G  O F  T H E  A L G O R I T H M" << std::endl << std::endl;    
    file << "  Number of different temperatures : " << config->nbt << std::endl;
    file << "  Temperature max : " << config->tmax << std::endl;
    file << "  Number of iterations " << config->nit << std::endl;
    file << "  Probability of allowing swapping : " << config->pee << std::endl;
    file << "  Number of components to modify in a MH iteration as a function of T : " << std::endl;
    for (int i=0;i<(int)(config->nc.size());i++)
      file << "    For T" << i << " : " << config->nc[i] << "  ";
    file << std::endl;
    file << "  Temperature ladder : ";
    for (int i=0;i<config->nbt;i++) 
      file << "T[" << i << "] = " << config->T[i] << "  ";
    file << std::endl;
    if (config->waveletParameterization) {
    
      file << std::endl << "W A V E L E T  C O N F I G" << std::endl << std::endl;
      file << "  Wavelet used : " << config->wavelet << std::endl;
      file << "  Number of DWT stages : " << config->ndwts << std::endl;
      for (int wavelet = 0; wavelet < config->nWavelets; wavelet++) { // Loop on the wavelets used
        file << "  Size of coeffsP[" << config->listOfWavelets[wavelet] << "] : " << config->coeffsP[wavelet].size() << std::endl;
        file << "  Size of flagP[" << config->listOfWavelets[wavelet] << "] : " << config->flagP[wavelet].size() << std::endl;
        file << "  Size of lengthP[" << config->listOfWavelets[wavelet] << "] : " << config->lengthP[wavelet].size() << std::endl;
        file << "  Size of coeffsS[" << config->listOfWavelets[wavelet] << "] : " << config->coeffsS[wavelet].size() << std::endl;
        file << "  Size of flagS[" << config->listOfWavelets[wavelet] << "] : " << config->flagS[wavelet].size() << std::endl;
        file << "  Size of lengthS[" << config->listOfWavelets[wavelet] << "] : " << config->lengthS[wavelet].size() << std::endl;
      }
      file << std::endl;
    }
    file.close();
  }
}

void copyDataFiles(const Configuration* config)
// Copy the data files used during the simulation in OUTPUT_FILES
{
  if (config->mpiConfig.rank == 0) { // In case of MPI just one process must manage files
    // Copy the configuration file
    std::string nameOfConfFile = NAME_OF_CONFIGURATION_FILE;
    std::string command="cp "+nameOfConfFile+" "+config->outputDir+nameOfConfFile;
    int status = system(command.c_str());
    if(status != 0) {
      std::cout << "Error while copying "+nameOfConfFile+" to ";
      std::cout << config->outputDir+nameOfConfFile+" (status = "<< status << ")" << std::endl;
      exit(1);
    }
    // Copy stations file
    command="cp "+config->filesDir+config->name_of_stations_file+" "+config->outputDir+config->name_of_stations_file;
    status = system(command.c_str());
    if(status != 0) {
      std::cout << "Error while copying "+config->filesDir+config->name_of_stations_file+" to ";
      std::cout << config->outputDir+config->name_of_stations_file+" (status = "<< status << ")" << std::endl;
      exit(1);
    }
    // Copy shots file
    command="cp "+config->filesDir+config->name_of_shots_file+" "+config->outputDir+config->name_of_shots_file;
    status = system(command.c_str());
    if(status != 0) {
      std::cout << "Error while copying "+config->filesDir+config->name_of_shots_file+" to ";
      std::cout << config->outputDir+config->name_of_shots_file+" (status = "<< status << ")" << std::endl;
      exit(1);
    }
    // Copy prior features file
    if (!config->buildPrior) {
      command="cp "+config->filesDir+config->name_of_prior_features_file+" "+config->outputDir+config->name_of_prior_features_file;
      status = system(command.c_str());
      if(status != 0) {
        std::cout << "Error while copying "+config->filesDir+config->name_of_prior_features_file+" to ";
        std::cout << config->outputDir+config->name_of_prior_features_file+" (status = "<< status << ")" << std::endl;
        exit(1);
      }
    }
    // Copy P first guess file
    command="cp "+config->filesDir+config->name_of_first_guess_P_file+" "+config->outputDir+config->name_of_first_guess_P_file;
    status = system(command.c_str());
    if(status != 0) {
      std::cout << "Error while copying "+config->filesDir+config->name_of_first_guess_P_file+" to ";
      std::cout << config->outputDir+config->name_of_first_guess_P_file+" (status = "<< status << ")" << std::endl;
      exit(1);
    }
    if (config->swaves) {
      // Copy S first guess file
      command="cp "+config->filesDir+config->name_of_first_guess_S_file+" "+config->outputDir+config->name_of_first_guess_S_file;
      status = system(command.c_str());
      if(status != 0) {
        std::cout << "Error while copying "+config->filesDir+config->name_of_first_guess_S_file+" to ";
        std::cout << config->outputDir+config->name_of_first_guess_S_file+" (status = "<< status << ")" << std::endl;
        exit(1);
      }
    }
    if (config->analyticalRun) {
      // Copy real P wave velocity file
      command="cp "+config->filesDir+config->name_of_real_profile_P+" "+config->outputDir+config->name_of_real_profile_P;
      status = system(command.c_str());
      if(status != 0) {
        std::cout << "Error while copying "+config->filesDir+config->name_of_real_profile_P+" to ";
        std::cout << config->outputDir+config->name_of_real_profile_P+" (status = "<< status << ")" << std::endl;
        exit(1);
      }
      if (config->swaves) {
        // Copy real S wave velocity file
        command="cp "+config->filesDir+config->name_of_real_profile_S+" "+config->outputDir+config->name_of_real_profile_S;
        status = system(command.c_str());
        if(status != 0) {
          std::cout << "Error while copying "+config->filesDir+config->name_of_real_profile_S+" to ";
          std::cout << config->outputDir+config->name_of_real_profile_S+" (status = "<< status << ")" << std::endl;
          exit(1);
        }
      }
    }
    else {
      // Copy times file
      command="cp "+config->filesDir+config->name_of_times_file+" "+config->outputDir+config->name_of_times_file;
      status = system(command.c_str());
      if(status != 0) {
        std::cout << "Error while copying "+config->filesDir+config->name_of_times_file+" to ";
        std::cout << config->outputDir+config->name_of_times_file+" (status = "<< status << ")" << std::endl;
        exit(1);
      }
    }
  }
}

void printEvolution(const int n, const Run* run, const Configuration* config)
// Print the iteration number n. TODO : add a % done and at when the run will finish.
{
  if (config->mpiConfig.rank == 0 && n%config->view == 0) {
    if (n==0)
      std::cout << std::endl << "********************** MAIN LOOP ***********************" << std::endl << std::endl;
    else {
      std::cout << "ITERATION NUMBER : " << n << " / " << config->nit << " (" << (double)n*100.0/((double)config->nit) << " %)"<< std::endl;
    }
  }
}

void writeFiles(Run* run, Configuration* config, int n)
// Write informations on files
{
  writeSCI(run,config);    // Write the Importance Weights matrix SCI in a file (sci.XXX.dat)
  writeStats(run,config);  // Write the actual statistics of the run in files (statsI.XXX.dat) 
  writeStatus(run,config); // Write the status of the chains in files (paramsChainI.XXX.dat)
  writeAverages(run,config);  // Write the average, variance and quantiles profiles on files
  if (n%config->iterationsBestProfiles == 0 && n > 1)
    writeBestProfiles(run,config); // Write the best profiles on files
  writeMinMaxProfiles(run,config); // Update min and max velocities explored, write it on file 
  if (config->computeResiduals && n%config->iterationsResiduals == 0 && n>1)
    computeResidualsForBestModel(run,config);
}

void writeStatus(Run* run,Configuration* config)
// If it does not exist, creates the data files (chainI.XXX.dat). Writes a new line on it.
// TODO : write c++ fashion
//    std::string conf(config->outputDir+"config."+config->code+".dat");  // We will obtain /home/abottero/config.321.dat for example
//    std::ofstream fe(conf.c_str()); // Create the file
{
  if (config->mpiConfig.rank == 0) {       // (In case of parallel implementation just one process has to create files)
    for(int i=0; i<(int)run->chains.size();i++) { // Loop on all the chains
      std::ostringstream ii; // Store i as a string
      ii << i;
      std::string file; // To store the name of the file
      file = config->outputDir+"chain"+ii.str()+"."+config->code+".dat" ;  // We will obtain /home/abottero/chain.321.dat for example
      FILE* fx=fopen(file.c_str(),"a"); // Create the file if it does not exist and open it
      for(int k=0; k<(int)run->chains[i]->states.back().params.size();k++) { // Loop on all the parameters
        fprintf(fx,"%6.4f ",run->chains[i]->states.back().params[k]); // Write each parameter
      }
      fprintf(fx,"%6.4f\n",run->chains[i]->states.back().E); // Write the energy at the end of the line
      fclose(fx); // Close the file
    }
  }
}

void writeSwap(Run* run,Configuration* config)
// If it does not exist, creates the data file (ll.XXX.dat). Writes a new line on it. It will store the swapping history
// TODO : write c++ fashion
{
  if (config->mpiConfig.rank == 0) {   // (In case of parallel implementation just one process has to create files)
    std::string file; // To store the name of the file
    file = config->outputDir+"ll."+config->code+".dat" ;  // We will obtain /home/abottero/ll.321.dat for example
    FILE* ll=fopen(file.c_str(),"a"); // Create the file if it does not exist and open it 
    fprintf(ll,"%d %d %d\n",run->swapHist.back().idxState,run->swapHist.back().idxChain,run->swapHist.back().sp);
    fclose(ll); // Close the file
  }
}

void writeSCI(Run* run,Configuration* config)
// Write the Importance Weights matrix SCI in sci.XXX.dat.
// TODO : write c++ fashion
{
  if (config->mpiConfig.rank == 0) {   // (In case of parallel implementation just one process has to create files)
    std::string file; // To store the name of the file
    file = config->outputDir+"sci."+config->code+".dat";  // We will obtain /home/abottero/sci.321.dat for example
    FILE* sci=fopen(file.c_str(),"a"); // Create the file if it does not exist and open it 
    for(unsigned int k=0;k<run->chains.size()-1;k++)
      fprintf(sci,"%Le ",run->SCI.back()[k]);
    fprintf(sci,"\n");
    fclose(sci); // Close the file
  }
}

void writeStats(Run* run,Configuration* config)
// Write the statistics of the chains (probability of swapping...) in statsI.XXX.dat.
// In the order : at rt od ps as rs acceptance probability swapping probability : see structure Chain
// TODO : write c++ fashion
{
  if (config->mpiConfig.rank == 0) {   // (In case of parallel implementation just one process has to create files)
    for(int i=0; i<(int)run->chains.size();i++)  // Loop on all the chains
    {
      std::ostringstream ii;   // Store i as a string
      ii << i;
      std::string file; // To store the name of the file
      file = config->outputDir+"stats"+ii.str()+"."+config->code+".dat"; // We will obtain /home/abottero/stats1.321.dat for example
      FILE* stats=fopen(file.c_str(),"a"); // Create the file if it does not exist and open it 
      fprintf(stats,"%d %d %d %d %d %d %f %f\n",run->chains[i]->at,run->chains[i]->rt,run->chains[i]->od,run->chains[i]->ps,run->chains[i]->as,run->chains[i]->rs,(double)run->chains[i]->at*100/(run->chains[i]->at+run->chains[i]->rt),(double)run->chains[i]->as*100/run->chains[i]->ps);
      fclose(stats); // Close the file
    }
  }
}

void printStatus(Run* run,Configuration* config,int view)
// Print the status of the run, display the "view" last states
{
  if (config->mpiConfig.rank == 0) {   // (In case of parallel implementation just one process has to print messages)
    std::cout << std::endl;
    std::cout << "****************STATUS******************" << std::endl;
    std::cout << "Number of chains : " << run->chains.size() << std::endl;
    std::cout << std::endl;
    for(int i=0;i<(int)run->chains.size();i++) { // Loop on the number of chains on the run
      std::cout << "-----------" << std::endl;
      std::cout << "_Chain["<< i <<"] : " << std::endl;
      std::cout << std::endl;
      std::cout << "  Length of the chain : "<< run->chains[i]->states.size() << std::endl;
      std::cout << "  Temperature : "<< run->chains[i]->T << std::endl;
      std::cout << "  Number of components to modify at each step : "<< run->chains[i]->nc << std::endl << std::endl;
      std::cout <<"Three last steps : "<< std::endl ;
      std::cout <<"  X    Y      E"<< std::endl;
      for(int j=(int)run->chains[i]->states.size()-view; j<(int)run->chains[i]->states.size();j++) { // Loop on the last "view" states
        for(int k=0; k<(int)run->chains[i]->states[0].params.size();k++) { // Loop on all the parameters
          printf("%5.3f ",run->chains[i]->states[j].params[k]);
        }
        printf("%5.3f\n",run->chains[i]->states[j].E);
      }
      std::cout << std::endl;
      std::cout << "Number of accepted MH transitions : " << run->chains[i]->at << std::endl;
      std::cout << "Number of rejected MH transitions : " << run->chains[i]->rt << std::endl;
      std::cout << "Acceptance probability : " << (double)run->chains[i]->at*100/(run->chains[i]->at+run->chains[i]->rt)<< " % " << std::endl;
      std::cout << "Number of -out of domain- : " << run->chains[i]->od << std::endl;
      std::cout << "Number of proposed swaps : " << run->chains[i]->ps << std::endl;
      std::cout << "Number of accepted swaps : " << run->chains[i]->as << std::endl;
      std::cout << "Number of rejected swaps : " << run->chains[i]->rs << std::endl;
      std::cout << "Probability of swapping : : " << (double)run->chains[i]->as*100/run->chains[i]->ps << " % " <<std::endl;
      std::cout << "-----------" << std::endl;
    }
  //  std::cout << "Last 5 lines of the matrix of importance weights :" << std::endl;
  //  for(unsigned int i=run->SCI.size()-5; i<run->SCI.size();i++) {
  //      for(unsigned int k=0;k<run->chains.size()-1;k++) {
  //          std::cout << run->SCI[i][k] << " ";
  //      }
  //      std::cout <<  std::endl;
  //  }
    std::cout << "****************END OF STATUS******************" << std::endl;
  }
}

void summary(Run* run, const Configuration* config)
// Display a summary of the run
{
  if (config->mpiConfig.rank == 0) {   // (In case of parallel implementation just one process has to print messages)
    std::cout << std::endl << "**** SUMMARY OF THE RUN: **** "<< std::endl;
    std::cout << std::endl;
    for(int i=0;i<(int)run->chains.size();i++) { // Loop on the chains
      std::cout << "_Chain ["<< i <<"] " << std::endl;
      std::cout << "Number of accepted MH transitions : " << run->chains[i]->at << std::endl;
      std::cout << "Number of rejected MH transitions : " << run->chains[i]->rt << std::endl;
      std::cout << "Acceptance probability : " << (double)run->chains[i]->at*100/(run->chains[i]->at+run->chains[i]->rt)<< " % " << std::endl;
      std::cout << "Number of -out of domain- : " << run->chains[i]->od << std::endl;
      std::cout << "Number of proposed swaps : " << run->chains[i]->ps << std::endl;
      std::cout << "Number of accepted swaps : " << run->chains[i]->as << std::endl;
      std::cout << "Number of rejected swaps : " << run->chains[i]->rs << std::endl;
      std::cout << "Probability of swapping : : " << (double)run->chains[i]->as*100/run->chains[i]->ps << " % " <<std::endl;
      std::cout << "-----------" << std::endl;
    }
    std::cout << "***************************** "<< std::endl << std::endl;
  }
}

void write_summary(const Run* run, const Configuration* config)
// Write a summary of the run at the end of the file config.XXX.dat
{
  if (config->mpiConfig.rank == 0) {   // (In case of parallel implementation just one process has to create files)
    std::string conf(config->outputDir+"config."+config->code+".dat");  // We will obtain /home/abottero/config.321.dat for example
    std::ofstream file(conf.c_str(), std::ios::app); // Open the file and put the cursor at the end
    if(file == NULL) {
      std::cout << "IMPOSSIBLE TO OPEN config." << config->code << ".dat" << std::endl;
      exit(0);
    }
    file << std::endl;
    file << "S U M M A R Y  O F  T H E  R U N"<< std::endl;
    file << std::endl;
    for(int i=0;i<(int)run->chains.size();i++) {
      file << "  _Chain ["<< i <<"] " << std::endl;
      file << "    Number of accepted MH transitions : " << run->chains[i]->at << std::endl;
      file << "    Number of rejected MH transitions : " << run->chains[i]->rt << std::endl;
      file << "    Acceptance probability : " << (double)run->chains[i]->at*100/(run->chains[i]->at+run->chains[i]->rt)<< " % " << std::endl;
      file << "    Number of -out of domain- : " << run->chains[i]->od << std::endl;
      file << "    Number of proposed swaps : " << run->chains[i]->ps << std::endl;
      file << "    Number of accepted swaps : " << run->chains[i]->as << std::endl;
      file << "    Number of rejected swaps : " << run->chains[i]->rs << std::endl;
      file << "    Probability of swapping : : " << (double)run->chains[i]->as*100/run->chains[i]->ps << " % " <<std::endl;
      file << std::endl << std::endl;
    }
    file << std::endl;
    time_t t=time(NULL);
    if (t == -1) { // Sometimes time() returns -1
      std::cout << "TIME PROBLEM"<< std::endl;
      exit(0);
    }
    char buffer[80];
    convertTime(buffer, t);
    file << "Final time : " << buffer << std::endl;
    file.close();
  }
}

void writeAverages(const Run* run, const Configuration* config)
// Write the average, variance and quantiles profiles on files
{
  if (config->mpiConfig.rank == 0) {   // (In case of parallel implementation just one process has to create files)
    std::string strFileGlobalAverage,strFileGlobalVar;
    strFileGlobalVar = config->outputDir+"globalVarP."+config->code+".dat";   // We will obtain /home/abottero/globalVarP.451.dat for example
    strFileGlobalAverage = config->outputDir+"globalAverageP."+config->code+".dat";   // We will obtain /home/abottero/globalAverageP.451.dat for example
    std::ofstream fileGlobalVar(strFileGlobalVar.c_str()),fileGlobalAverage(strFileGlobalAverage.c_str());
    for (int iz=0;iz<config->data.nzFilt-1;iz++) { // Loop on the values
      fileGlobalAverage << config->data.zFiltp[iz] << " " << run->averageP[iz] << std::endl;
      fileGlobalVar << config->data.zFiltp[iz] << " " << run->varP[iz] << std::endl;
    }
    fileGlobalAverage.close(); 
    fileGlobalVar.close(); 
    if(config->swaves) {
      strFileGlobalVar = config->outputDir+"globalVarS."+config->code+".dat";   // We will obtain /home/abottero/globalVarS.451.dat for example
      strFileGlobalAverage = config->outputDir+"globalAverageS."+config->code+".dat";   // We will obtain /home/abottero/globalAverageS.451.dat for example
      std::ofstream fileGlobalVar(strFileGlobalVar.c_str()),fileGlobalAverage(strFileGlobalAverage.c_str());
      for (int iz=0;iz<config->data.nzFilt-1;iz++) { // Loop on the values
        fileGlobalVar << config->data.zFiltp[iz] << " " << run->varS[iz] << std::endl;
        fileGlobalAverage << config->data.zFiltp[iz] << " " << run->averageS[iz] << std::endl;
      }
      fileGlobalAverage.close(); 
      fileGlobalVar.close(); 
    }
    for(int i=0; i<(int)run->chains.size();i++) { // Loop on all the chains
      std::ostringstream ii;   // Store i as a string
      ii << i;
      std::string strFileAverage,strFileVar, strFileQinf, strFileQsup; // To store the name of the files
      strFileAverage = config->outputDir+"averageP"+ii.str()+"."+config->code+".dat"; // We will obtain /home/abottero/averageP0.451.dat for example
      strFileVar = config->outputDir+"varP"+ii.str()+"."+config->code+".dat";   // We will obtain /home/abottero/varP0.451.dat for example
      strFileQinf = config->outputDir+"qInfP"+ii.str()+"."+config->code+".dat"; // We will obtain /home/abottero/qInfP0.451.dat for example
      strFileQsup = config->outputDir+"qSupP"+ii.str()+"."+config->code+".dat"; // We will obtain /home/abottero/qSupP0.451.dat for example
      std::ofstream fileAverage(strFileAverage.c_str()),fileVar(strFileVar.c_str()),fileQinf(strFileQinf.c_str()),fileQsup(strFileQsup.c_str()); 
      for (int iz=0;iz<config->data.nzFilt-1;iz++) { // Loop on the values
        fileAverage << config->data.zFiltp[iz] << " " << run->chains[i]->averageP[iz] << std::endl;
        fileVar << config->data.zFiltp[iz] << " " << run->chains[i]->varP[iz]/((double)((int)run->chains[i]->profilesP[iz].size()-1)) << std::endl;
        fileQinf << config->data.zFiltp[iz] << " " << run->chains[i]->qInfP[iz] << std::endl;
        fileQsup << config->data.zFiltp[iz] << " " << run->chains[i]->qSupP[iz] << std::endl;
      }
      fileAverage.close(); 
      fileVar.close(); 
      fileQinf.close(); 
      fileQsup.close(); 
      
      if(config->swaves) {
        strFileAverage = config->outputDir+"averageS"+ii.str()+"."+config->code+".dat"; // We will obtain /home/abottero/averageS0.451.dat for example
        strFileVar = config->outputDir+"varS"+ii.str()+"."+config->code+".dat";   // We will obtain /home/abottero/varS0.451.dat for example
        strFileQinf = config->outputDir+"qInfS"+ii.str()+"."+config->code+".dat"; // We will obtain /home/abottero/qInfS0.451.dat for example
        strFileQsup = config->outputDir+"qSupS"+ii.str()+"."+config->code+".dat"; // We will obtain /home/abottero/qSupS0.451.dat for example
        std::ofstream fileAverageS(strFileAverage.c_str()),fileVarS(strFileVar.c_str()),fileQinfS(strFileQinf.c_str()),fileQsupS(strFileQsup.c_str()); 
        for (int iz=0;iz<config->data.nzFilt-1;iz++) { // Loop on the values
          fileAverageS << config->data.zFiltp[iz] << " " << run->chains[i]->averageS[iz] << std::endl;
          fileVarS << config->data.zFiltp[iz] << " " << run->chains[i]->varS[iz]/((double)((int)run->chains[i]->profilesS[iz].size()-1))  << std::endl;
          fileQinfS << config->data.zFiltp[iz] << " " << run->chains[i]->qInfS[iz] << std::endl;
          fileQsupS << config->data.zFiltp[iz] << " " << run->chains[i]->qSupS[iz] << std::endl;
        }
        fileAverageS.close(); 
        fileVarS.close(); 
        fileQinfS.close(); 
        fileQsupS.close();       
      }
    }
  }
}

void writeBestProfiles(Run* run, const Configuration* config)
// Write the best profiles on files 
// TODO : This could be far more optimized!
{
  std::vector<double>::const_iterator EmaxIt = std::max_element(run->bestE.begin(), run->bestE.end());
  std::vector<double> temp_profileP, temp_profileS;
  double Emax = *EmaxIt;
  int EmaxIdx = (int)(EmaxIt - run->bestE.begin());
  for(int i=0; i<(int)run->chains.size();i++) { // Loop on all the chains
  // std::cout << std::endl << "Chain " << i << ":" << std::endl;
    for(int j=0; j<(int)run->chains[i]->states.size();j++) { // Loop on all the states
    //  std::cout << run->chains[i]->states[j].E*run->chains[i]->T << " ";
      bool alreadyContained = false;
      for (unsigned int k=0; k<run->bestE.size();k++) {
        if (fabs(run->bestE[k]-run->chains[i]->states[j].E*run->chains[i]->T)<TINYVAL)
          alreadyContained = true;
      }
      if(!alreadyContained) { // run->bestE does not contains the state's energy
        if((int)run->bestE.size() < config->nBestProfiles) { // If size of run->bestE is smaller that config->nBestProfile
          run->bestE.push_back(run->chains[i]->states[j].E*run->chains[i]->T);
          run->idxE.push_back(j);
          run->chainBestE.push_back(i);
          EmaxIt = std::min_element(run->bestE.begin(), run->bestE.end());
          Emax = *EmaxIt; 
          EmaxIdx = (int)(EmaxIt - run->bestE.begin());
        }
        else if (run->chains[i]->states[j].E*run->chains[i]->T < Emax) { // If this state's energy is greater than the ones stored
          run->bestE[EmaxIdx] = run->chains[i]->states[j].E*run->chains[i]->T;
          run->idxE[EmaxIdx]=j;
          run->chainBestE[EmaxIdx]=i;
          EmaxIt = std::max_element(run->bestE.begin(), run->bestE.end());
          Emax = *EmaxIt; 
          EmaxIdx = (int)(EmaxIt - run->bestE.begin());
        }
      }
    }
  //  std::cout << std::endl;
  }
 // std::cout << "Best energies :" << std::endl;
  int status;
  if (config->mpiConfig.rank == 0) {
    status = system("rm -f best*.dat"); // To remove the old files (this is not optimum at all) TODO
    if(status != 0) {
      std::cout << "Error while deleting best profiles files" << std::endl;
      exit(1);
    }
  }
  for (unsigned int k=0; k<run->bestE.size();k++) { // Loop on the best models
  //  std::cout << run->bestE[k] << " (chain : " << run->chainBestE[k] << " idx : " << run->idxE[k] << ")  ";
    temp_profileP.clear();
    temp_profileS.clear();
    for (int iz=0;iz<config->data.nzFilt-1;iz++) {
      temp_profileP.push_back(run->chains[run->chainBestE[k]]->profilesP[iz][run->idxE[k]]);
      if(config->swaves)
        temp_profileS.push_back(run->chains[run->chainBestE[k]]->profilesS[iz][run->idxE[k]]);
    }
    if (config->mpiConfig.rank == 0) {   // (In case of parallel implementation just one process has to create files)
      std::ostringstream iiChain, iiIdx, iiE;   // Store as strings
      iiChain << run->chainBestE[k];
      iiIdx << run->idxE[k];
      iiE << run->bestE[k];
      std::string nameP = config->outputDir+"bestPprofile.chain"+iiChain.str()+".idx"+iiIdx.str()+".E"+iiE.str()+"."+config->code+".dat"; 
      write_two_columns_file(&config->data.zFiltp,&temp_profileP, nameP);
      if(config->swaves) {
        std::string nameS = config->outputDir+"bestSprofile.chain"+iiChain.str()+".idx"+iiIdx.str()+".E"+iiE.str()+"."+config->code+".dat"; 
        write_two_columns_file(&config->data.zFiltp,&temp_profileS, nameS);
      }
    }
  }
}

void computeResidualsForBestModel(Run* run, Configuration* config)
// Recompute the eikonal (in serial) for the best model and save the residuals on a file
{
  if (config->mpiConfig.rank == 0) 
    std::cout << "Recomputing forward problem (in serial) for best model..." << std::endl;    
  int temp_nBestProfiles=config->nBestProfiles;
  config->nBestProfiles=1; // If the user have set 0 the following would not work
  writeBestProfiles(run, config); // Store the properties of the best profile in run->bestE[k],run->idxE[k], run->chainBestE[k]
  // write the best profile on a file
  std::vector<double> bestProfileP, bestProfileS;
  std::vector<double>::const_iterator EminIt = std::min_element(run->bestE.begin(), run->bestE.end());
  int EminIdx = (int)(EminIt - run->bestE.begin()); // Index of the minimum energy
  // double Emin = *EminIt;
  for (int iz=0;iz<config->data.nzFilt-1;iz++) {
    bestProfileP.push_back(run->chains[run->chainBestE[EminIdx]]->profilesP[iz][run->idxE[EminIdx]]);
    if(config->swaves)
      bestProfileS.push_back(run->chains[run->chainBestE[EminIdx]]->profilesS[iz][run->idxE[EminIdx]]);
  }
  VelocityModel velModel;
  velModel.nx=config->data.nx; velModel.ny=config->data.ny; velModel.nz=config->data.nzFilt;
  velModel.dx=config->data.dx; velModel.dy=config->data.dy; velModel.dz=config->data.dzFilt;
  velModel.xmin = config->data.xminGrid; velModel.ymin = config->data.yminGrid; velModel.zmin = config->data.zminGrid;
  velModel.velP= new tab3d<double>(velModel.nz-1,velModel.nx-1,velModel.ny-1,-1.0); // Create a (nz-1,nx-1,ny-1) mesh and initialize every cell at -1 
 // Copy the file content into the velocity model (extend the 1D profile to obtain a 3D profile).
  meshing(&bestProfileP,&velModel,false); // Extend this profile on the whole mesh 
  if (config->swaves) {
    velModel.velS= new tab3d<double>(velModel.nz-1,velModel.nx-1,velModel.ny-1,-1.0); // Create a (nz-1,nx-1,ny-1) mesh and initialize every cell at -1 
    // Copy the file content into the velocity model (extend the 1D profile to obtain a 3D profile).
    meshing(&bestProfileS,&velModel,true); // Extend this profile on the whole mesh
  } 
  tab3d<double> tt3dP(config->data.nzFilt,config->data.nx,config->data.ny,-1.0); // Will contain the P waves arrival times. Initialize every cell at -1.0 
  tab3d<double> tt3dS(config->data.nzFilt,config->data.nx,config->data.ny,-1.0); // TODO : this is declared even if swaves == false -> change that
//  double E=0.0;
  ArrivalTimes arrivalTimes;
  int numberOfShots = (int)config->data.coordShots.size();
  int numberOfStations = (int)config->data.coordStations.size();
  std::vector<int> indexP, indexS;
 // double sumP=0.0, sumS=0.0;
  double t0=0.0;   // Origin time of the shot (if config->recalculateT0 == 1 we will recalculate it)
  int numberOfEikonalToCompute = 0;
  if(config->swaves) 
    numberOfEikonalToCompute = 2*numberOfShots;
  else
    numberOfEikonalToCompute = numberOfShots;
   
  for(int i=0;i<numberOfEikonalToCompute;i++) { // Loop on the shots
    if (i<numberOfShots) { // P waves
      if(config->verbose2 && config->mpiConfig.rank == 0)
        std::cout << "     Eikonal computing for best P wave velocity profile, shot number " << i+1 << " on " << numberOfShots << " ..."<< std::endl;
      indexP.push_back(i); 
      eikonal3d(&tt3dP,&velModel,config,i,false); 
      //Calculate the P waves travel times everywhere on the mesh (put them on tt3dP) for the shot number i
      for(int j=0;j<numberOfStations;j++) 
        arrivalTimes.timesP.push_back(getTime(&tt3dP,config->data.coordStations[j],&velModel));

      if (config->recalculateT0) {
        int nUsedP=0;  // Number of travel times used for P waves
        for(int k=0;k<(int)arrivalTimes.timesP.size();k++) { // Calculation of t0
          int idxP=indexP[(int)floor(k/numberOfStations)]*numberOfStations+k-numberOfStations*(int)floor(k/numberOfStations); // index of the arrival time
          if (config->data.times.timesP[idxP] > 0.0) {
            t0+=config->data.times.timesP[idxP]-arrivalTimes.timesP[k];
            nUsedP++;
          }
        }
        t0=t0/nUsedP;
      }
    } 
    else { // S waves
      if(config->verbose2 && config->mpiConfig.rank == 0)
        std::cout << "     Eikonal computing for best S wave velocity profile, shot number " << i+1-numberOfShots << " on " << numberOfShots << " ..."<< std::endl;
      indexS.push_back(i-numberOfShots);
      eikonal3d(&tt3dS,&velModel,config,i-numberOfShots,true); 
      // Calculate the S waves travel times everywhere on the mesh (put them on tt3dS) for the shot number i
      for(int j=0;j<numberOfStations;j++) 
        arrivalTimes.timesS.push_back(getTime(&tt3dS,config->data.coordStations[j],&velModel));
    }  
  }

/*
  for(int k=0;k<(int)arrivalTimes.timesP.size();k++) { // Loop on the P wave travel times calculated
    int idxP=indexP[(int)floor(k/numberOfStations)]*numberOfStations+k-numberOfStations*(int)floor(k/numberOfStations); // index of the arrival time corresponding to the one calculated. 
    double diffP=0.0;
    if (config->data.times.timesP[idxP] > 0.0)
      diffP= config->data.times.timesP[idxP]-(t0+arrivalTimes.timesP[k]);
 //   std::cout << " config->data.times.timesP[" << idxP << "]  : " << config->data.times.timesP[idxP] << "  arrivalTimes.timesP[" << k << "] : " << arrivalTimes.timesP[k] << "  Diff : " << diffP << std::endl; 
 //   sumP+=pow(diffP/config->data.sigmaP,2.0)/2.0; // sum of the squares of the differences between the P wave first arrival times and the data
  }
  if(config->swaves) { 
    for(int k=0;k<(int)arrivalTimes.timesS.size();k++) { // Loop on the S wave travel times calculated
      int idxS=indexS[(int)floor(k/numberOfStations)]*numberOfStations+k-numberOfStations*(int)floor(k/numberOfStations); // index of the arrival time corresponding to the one calculated
      double diffS=0.0;
      if (config->data.times.timesS[idxS] > 0.0)
        diffS= config->data.times.timesS[idxS]-(t0+arrivalTimes.timesS[k]);
    //    std::cout << " config->data.times.timesS[" << idxS << "]  : " << config->data.times.timesS[idxS] << "  arrivalTimes.timesS[" << k << "] : " << arrivalTimes.timesS[k] << "  Diff : " << diffS << std::endl;   
 //     sumS+=pow(diffS/config->data.sigmaS,2.0)/2.0; // sum of the squares of the differences between the P wave first arrival times
    }
  }
//  E=(sumP+sumS+config->data.Ep)/(run->chains[run->chainBestE[EminIdx]]->T);
//  std::ostringstream ii;   // Store E as a string
//  ii << E;
*/
  if(config->mpiConfig.rank == 0) {
    if(config->swaves)
      write_two_columns_file(&arrivalTimes.timesP,&arrivalTimes.timesS, config->outputDir+"bestModelTimes."+config->code+".dat");
    else
      write_one_column_file(&arrivalTimes.timesP, config->outputDir+"bestModelTimes."+config->code+".dat");
  }
  delete velModel.velP;
  if(config->swaves) 
    delete velModel.velS;
  config->nBestProfiles=temp_nBestProfiles;
  if (config->mpiConfig.rank == 0) 
    std::cout << "Done !" << std::endl;
}

void writeMinMaxProfiles(Run* run, const Configuration* config)
// Update min and max velocities explored by the run, write them on file
{

  if(config->verbose2 && config->mpiConfig.rank == 0)
    std::cout << "Updating max and min profiles investigated by the run... " << std::endl;
  int it=run->chains[0]->profilesP[0].size()-1;
  
  //************** Update MIN and MAX profiles *****************//
  double lastPvel = 0.0, lastSvel=0.0;
  
  for(int i=0;i<config->nbt;i++) { // Loop on the chains
    for(int iz=0;iz<config->data.nzFilt-1;iz++) { // Loop on the depths
      lastPvel = run->chains[i]->profilesP[iz].back();
      if (lastPvel > run->chains[i]->maxP[iz] || it < TRESH)
        run->chains[i]->maxP[iz] = lastPvel;
      if (lastPvel < run->chains[i]->minP[iz] || it < TRESH)
        run->chains[i]->minP[iz] = lastPvel;
      if (config->swaves) {
        lastSvel = run->chains[i]->profilesS[iz].back();
        if (lastSvel > run->chains[i]->maxS[iz] || it < TRESH)
          run->chains[i]->maxS[iz] = lastSvel;
        if (lastSvel < run->chains[i]->minS[iz] || it < TRESH)
          run->chains[i]->minS[iz] = lastSvel;
      }
    }
    if (config->mpiConfig.rank == 0) {   // (In case of parallel implementation just one process has to create files)
      // ************** write them on files ***************** //
      std::ostringstream ii;   // Store i as a string
      ii << i;
      write_two_columns_file(&config->data.zFiltp,&run->chains[i]->maxP, config->outputDir+"maxP."+ii.str()+"."+config->code+".dat");
      write_two_columns_file(&config->data.zFiltp,&run->chains[i]->minP, config->outputDir+"minP."+ii.str()+"."+config->code+".dat");
      if(config->swaves) {
        write_two_columns_file(&config->data.zFiltp,&run->chains[i]->maxS, config->outputDir+"maxS."+ii.str()+"."+config->code+".dat");
        write_two_columns_file(&config->data.zFiltp,&run->chains[i]->minS, config->outputDir+"minS."+ii.str()+"."+config->code+".dat");
      }
    }
  }
  std::vector<double> minPforThisDepth,maxPforThisDepth,minSforThisDepth,maxSforThisDepth;
  for(int i=0;i<config->nbt;i++) { // Loop on the chains
    minPforThisDepth.push_back(0.0);
    maxPforThisDepth.push_back(0.0);
    minSforThisDepth.push_back(0.0);
    maxSforThisDepth.push_back(0.0); 
  }
  for(int iz=0;iz<config->data.nzFilt-1;iz++) { // Loop on the depths
    for(int i=0;i<config->nbt;i++) { // Loop on the chains
      minPforThisDepth[i]=run->chains[i]->minP[iz];
      maxPforThisDepth[i]=run->chains[i]->maxP[iz];
      minSforThisDepth[i]=run->chains[i]->minS[iz];
      maxSforThisDepth[i]=run->chains[i]->maxS[iz];
    }
    if(it > TRESH) {
      run->minP[iz] = *std::min_element(minPforThisDepth.begin(), minPforThisDepth.end());
      run->maxP[iz] = *std::max_element(maxPforThisDepth.begin(), maxPforThisDepth.end());
      run->minS[iz] = *std::min_element(minSforThisDepth.begin(), minSforThisDepth.end());
      run->maxS[iz] = *std::max_element(maxSforThisDepth.begin(), maxSforThisDepth.end());  
    }
    else {
      run->minP[iz] = run->chains[0]->minP[iz];
      run->maxP[iz] = run->chains[0]->maxP[iz];
      run->minS[iz] = run->chains[0]->minS[iz];
      run->maxS[iz] = run->chains[0]->maxS[iz];  
    }
  }
  if (config->mpiConfig.rank == 0) {   // (In case of parallel implementation just one process has to create files)
    write_two_columns_file(&config->data.zFiltp,&run->maxP, config->outputDir+"maxP."+config->code+".dat");
    write_two_columns_file(&config->data.zFiltp,&run->minP, config->outputDir+"minP."+config->code+".dat");
    if(config->swaves) {
      write_two_columns_file(&config->data.zFiltp,&run->maxS, config->outputDir+"maxS."+config->code+".dat");
      write_two_columns_file(&config->data.zFiltp,&run->minS, config->outputDir+"minS."+config->code+".dat");
    }
  }
  if(config->verbose2 && config->mpiConfig.rank == 0)
    std::cout << "Done!" << std::endl << std::endl;
}

void write_one_column_file(const std::vector<double>* column, const std::string name_of_file)
// Write a vector into a column file
{
  std::ofstream output_file(name_of_file.c_str());
  if (output_file.is_open()) {
    for (unsigned int i = 0; i <  (*column).size(); i++) // Loop on the coeffs
      output_file << std::setprecision(12) << (*column)[i] << std::endl;
  }
  else {
    std::cout << "Unable to open file "+name_of_file << std::endl;
    exit(0);
  }
}

void write_two_columns_file(const std::vector<double>* column1, const std::vector<double>* column2, const std::string name_of_file)
// Write a two columns file from two vectors of the same size
{
  std::ofstream output_file(name_of_file.c_str());
  if (output_file.is_open()) {
    if ((*column1).size() != (*column2).size())
      std::cout << "Impossible to create the file : "+name_of_file+" -> the two columns don't have the same size" << std::endl;
    for (unsigned int i = 0; i <  (*column1).size(); i++) // Loop on the coeffs
      output_file << std::setprecision(12) << (*column1)[i] << " " << (*column2)[i] << std::endl;
  }
  else {
    std::cout << "Unable to open file "+name_of_file << std::endl;
    exit(0);
  }
}

void write_three_columns_file(const std::vector<int>* column1, const std::vector<double>* column2, const std::vector<double>* column3, const std::string name_of_file)
// Write a three columns file from three vectors of the same size -> just used for prior features
{
  std::ofstream output_file(name_of_file.c_str());
  if (output_file.is_open()) {
    if (((*column1).size() != (*column2).size()) || ((*column1).size() != (*column3).size()))
      std::cout << "Impossible to create the file : "+name_of_file+" -> the three columns don't have the same size" << std::endl;
    for (unsigned int i = 0; i <  (*column1).size(); i++) // Loop on the coeffs
      output_file << (*column1)[i] << " " << std::setprecision(12) << (*column2)[i] << " " << (*column3)[i] << std::endl;
  }
  else {
    std::cout << "Unable to open file "+name_of_file << std::endl;
    exit(0);
  }
}

