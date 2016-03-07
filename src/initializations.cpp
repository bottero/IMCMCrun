#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <iomanip>      // std::setprecision
#include <algorithm>    // std::min_element, std::max_element
#include <numeric>      // accumulate
#include <stdlib.h>     // rand
#include "rvgs.h"
#include "rngs.h"
#include "defines.h"
#include "configFileParser.h"
#include "initializations.h"
#include "functions.h"
#include "structures.h"
#include "generalFunctions.h"
#include "filesAndControl.h"
#include "wavelet2s.h"

Configuration parameters_setting(int argc, char *argv[])
// Initialize the configuration structure containing the config of the run, read the configuration file
{
  Configuration config; // Create a structure config to store the configuration
  #ifdef PAR // Preprocessor instruction : if PAR is defined <-> if MPI mode (see in define.h)
    config.mpiConfig.rc=MPI_Init( &argc, &argv); // Starting MPI...
    if (config.mpiConfig.rc != MPI_SUCCESS) {    // Check if MPI has been lauch correctly
      std::cout << std::endl << "Error starting MPI program. Terminating" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, config.mpiConfig.rc); // Abort MPI
      exit(0); // Exit program
    }
    // MPI_Status status;
    MPI_Comm_rank( MPI_COMM_WORLD, &config.mpiConfig.rank); // Store rank of the process  
    MPI_Comm_size(MPI_COMM_WORLD, &config.mpiConfig.nb_process); // Store the number of processes
    MPI_Get_processor_name(config.mpiConfig.hostname, &config.mpiConfig.len); // Store processor name
    MPI_Barrier(MPI_COMM_WORLD); // Every process wait for the others here
  #else // If we are not in MPI :
    config.mpiConfig.rank=0; // rank = 0
    config.mpiConfig.nb_process=1; // and just one process working
  #endif
  config.startTime=time(NULL);
  if (config.startTime == -1) { // Sometimes time() returns -1
    std::cout << "TIME PROBLEM"<< std::endl;
    exit(0);
  }
  std::ostringstream startTime;   // Store a part of the starting time to distinguish files from different runs
  startTime << (int)config.startTime%1000; // -> gives 3 digit that will different each time we run the program
  config.code = startTime.str(); // Store these 3 digits on a string : config.code
  char buffer[80];
  convertTime(buffer, config.startTime);
  if (config.mpiConfig.rank == 0) {
    std::cout << "========================================================" << std::endl;
    std::cout << "                 IMPORTANCE RESAMPLING MCMC             " << std::endl;
    std::cout << "========================================================" << std::endl;
    std::cout << "This is the run number "<< config.code << std::endl;
    std::cout << "Starting time : "<< buffer << std::endl<< std::endl;
    #ifdef PAR  
      std::cout << "                  === MPI MODE ===" << std::endl;
      std::cout << config.mpiConfig.nb_process << " processes running on : " << config.mpiConfig.hostname << std::endl;
    #endif
    std::cout << "******************** INITIALIZATION ********************" << std::endl << std::endl;
  }
  std::string nameOfConfFile = NAME_OF_CONFIGURATION_FILE;
  readConfFile(nameOfConfFile,&config); // Read the informations contained in the configuration file
  if (!config.useAllWavelets or !config.waveletParameterization) {// If we are just using one wavelet...
    config.nWavelets = 1;
    config.listOfWavelets[0] = config.wavelet; // We store it at the first position of config.listOfWavelets
  }
  else
    config.nWavelets = sizeof(config.listOfWavelets)/sizeof(config.listOfWavelets[0]); // == NUMBER_OF_WAVELETS
  std::vector<double> emptyVectorOfDouble;
  std::vector<int> emptyVectorOfInt;
  for (int wavelet = 0; wavelet < config.nWavelets; wavelet++) { // Loop on the wavelets used
    config.coeffsP.push_back(emptyVectorOfDouble);
    config.coeffsS.push_back(emptyVectorOfDouble);
    config.flagP.push_back(emptyVectorOfDouble);
    config.flagS.push_back(emptyVectorOfDouble);
    config.lengthP.push_back(emptyVectorOfInt);
    config.lengthS.push_back(emptyVectorOfInt);
    config.data.filtFirstGuessP.push_back(emptyVectorOfDouble);
    config.data.filtFirstGuessS.push_back(emptyVectorOfDouble);
    config.data.minParameters.push_back(emptyVectorOfDouble);
    config.data.maxParameters.push_back(emptyVectorOfDouble);
    config.data.indexParameters.push_back(emptyVectorOfInt);
  }
  config.data.nz=-9999; // Will be updated in loadFirstGueses
  config.data.dz=-9999; // Will be updated in loadFirstGueses
  config.outputDir="./OUTPUT_FILES/"+config.code+"/";
  std::string mkdir_command="mkdir -p "+config.outputDir;
  int status;
  if (config.mpiConfig.rank == 0) {
    status = system(mkdir_command.c_str()); // TODO : This is a little bit dirty but there is no simple equivalent to "mkdir -p"
    if(status != 0) {
      std::cout << "Error while creating "+config.outputDir+" (status = "<< status << ")" << std::endl;
      exit(1);
    }
  }
  int seed=config.defaultSeed;
  if (config.useDefaultSeed != 1 && config.mpiConfig.rank == 0)
    seed=time(NULL);
  #ifdef PAR  
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&seed,1,MPI_DOUBLE,0,MPI_COMM_WORLD); 
    // The process 0 diffuses the seed to every other process
  #endif
  config.seed=seed;
  PlantSeeds(config.seed); // Seed initialization
  srand(config.seed);      // Seed initialization (useless?)
  build_temperature_ladders(&config); // Build the temperature ladders (T[0]=1)
  if (config.verbose1 && config.mpiConfig.rank == 0)
    std::cout << "Configuration structure initialized" << std::endl;
  return config;
}

void readConfFile(std::string nameOfConfFile, Configuration* config)
// Read the informations contained in the configuration file
{
  configFileParser::data myconfigdata;
  std::ifstream f(nameOfConfFile.c_str()); // Open the configuration file
  f >> myconfigdata; // store the content of the file in a configFileParser::data
  f.close(); // Close the file
  //std::cout << "This is the configuration file cleaned :" << std::endl;
  //std::cout << myconfigdata;
  //std::cout << std::endl;
  configFileParser::data::const_iterator iter;
  std::string tempNC,tempString;
  for (iter = myconfigdata.begin(); iter != myconfigdata.end(); iter++) {
    if (iter->first == "DATA_DIRECTORY")
      config->filesDir = iter->second;
    if (iter->first == "NPU")
      config->npu = atoi(iter->second.c_str());
    if (iter->first == "A_PRIOR")
      config->A = atof(iter->second.c_str());
    if (iter->first == "N_PRIOR_PROFILES")
      config->nPriorProfiles = atoi(iter->second.c_str());
    if (iter->first == "SIGMAP")
      config->data.sigmaP = atof(iter->second.c_str());
    if (iter->first == "SIGMAS")
      config->data.sigmaS = atof(iter->second.c_str());
    if (iter->first == "DI")
      config->di = atof(iter->second.c_str());
    if (iter->first == "DF")
      config->df = atof(iter->second.c_str());
    if (iter->first == "QP")
      config->qp = atof(iter->second.c_str());
    if (iter->first == "WAVELET_PARAMETERIZATION")
      config->waveletParameterization = atoi(iter->second.c_str());
    if (iter->first == "USE_ALL_WAVELETS")
      config->useAllWavelets = atoi(iter->second.c_str());
    if (iter->first == "KEEP_FIRST_VALUES")
      config->keep_first_values = atoi(iter->second.c_str());
    if (iter->first == "NAME_OF_REAL_PROFILE_FILE_P")
      config->name_of_real_profile_P = iter->second;
    if (iter->first == "NAME_OF_REAL_PROFILE_FILE_S")
      config->name_of_real_profile_S = iter->second;
    if (iter->first == "NAME_OF_FIRST_GUESS_P_FILE")
      config->name_of_first_guess_P_file = iter->second;
    if (iter->first == "NAME_OF_FIRST_GUESS_S_FILE")
      config->name_of_first_guess_S_file = iter->second;
    if (iter->first == "NAME_OF_TIMES_FILE")
      config->name_of_times_file = iter->second;
    if (iter->first == "NAME_OF_STATIONS_FILE")
      config->name_of_stations_file = iter->second; 
    if (iter->first == "NAME_OF_SHOTS_FILE")
      config->name_of_shots_file = iter->second;
    if (iter->first == "NAME_OF_PRIOR_FEATURES_FILE")
      config->name_of_prior_features_file = iter->second;       
    if (iter->first == "NSWEEPS")
      config->nSweeps = atoi(iter->second.c_str());
    if (iter->first == "EPSIN")
      config->epsin = atof(iter->second.c_str());
    if (iter->first == "NBT")
      config->nbt = atoi(iter->second.c_str());
    if (iter->first == "TMAX")
      config->tmax = atof(iter->second.c_str()); 
    if (iter->first == "NIT")
      config->nit = atoi(iter->second.c_str()); 
    if (iter->first == "PEE")
      config->pee = atof(iter->second.c_str()); 
    if (iter->first == "NC")
      tempNC = iter->second; 
    if (iter->first == "WAVELET")
      config->wavelet = iter->second;
    if (iter->first == "NDWTS")
      config->ndwts = atoi(iter->second.c_str());
    if (iter->first == "SWAVES")
      config->swaves = atoi(iter->second.c_str()); 
    if (iter->first == "VERBOSE1")
      config->verbose1 = atoi(iter->second.c_str());
    if (iter->first == "VERBOSE2")
      config->verbose2 = atoi(iter->second.c_str());
    if (iter->first == "VIEW")
      config->view = atoi(iter->second.c_str()); 
    if (iter->first == "COORD_TOL")
      config->coordTol = atof(iter->second.c_str());
    if (iter->first == "BUILD_PRIOR")
      config->buildPrior = atoi(iter->second.c_str());
    if (iter->first == "FIND_OPTIMUM_GRID")
      config->findOptimumGrid = atoi(iter->second.c_str());
    if (iter->first == "SHOT_NUMBER_REF")
      config->shotNumberRef = atoi(iter->second.c_str());
    if (iter->first == "NXREF")
      config->nxref = atoi(iter->second.c_str());
    if (iter->first == "NYREF")
      config->nyref = atoi(iter->second.c_str());
    if (iter->first == "NX_DEFAULT")
      config->nxDefault = atoi(iter->second.c_str());
    if (iter->first == "NY_DEFAULT")
      config->nyDefault = atoi(iter->second.c_str());
    if (iter->first == "NZFILT_DEFAULT")
      config->nzfiltDefault = atoi(iter->second.c_str());
    if (iter->first == "ANALYTICAL_RUN")
      config->analyticalRun = atoi(iter->second.c_str());
    if (iter->first == "RECALCULATE_T0")
      config->recalculateT0 = atoi(iter->second.c_str());
    if (iter->first == "USE_DEFAULT_SEED")
      config->useDefaultSeed = atoi(iter->second.c_str());
    if (iter->first == "DEFAULT_SEED")
      config->defaultSeed = atoi(iter->second.c_str());
    if (iter->first == "TEST")
      config->test = atoi(iter->second.c_str());
    if (iter->first == "MIN_E_TEST")
      config->minEtest = atof(iter->second.c_str());
    if (iter->first == "MAX_E_TEST")
      config->maxEtest = atof(iter->second.c_str());
    if (iter->first == "N_BEST_PROFILES")
      config->nBestProfiles = atoi(iter->second.c_str());
    if (iter->first == "COMPUTE_RESIDUALS")
      config->computeResiduals = atoi(iter->second.c_str());
    if (iter->first == "ITERATIONS_RESIDUALS")
      config->iterationsResiduals = atoi(iter->second.c_str());
    if (iter->first == "ITERATIONS_BEST_PROFILES")
      config->iterationsBestProfiles = atoi(iter->second.c_str());
    if (iter->first == "ONLY_CALCULATE_TIMES_FOR_FIRST_GUESS")
      config->calculateTimesForFirstGuess = atoi(iter->second.c_str());
    if (iter->first == "RESAMPLE")
      config->resample = atoi(iter->second.c_str());
    if (iter->first == "NXVEC") {
      tempString = trim(iter->second);
      int nValGiven = std::count(tempString.begin(), tempString.end(), ',') + 1;
      std::replace(tempString.begin(), tempString.end(), ',', ' ' ); // replace the comas with spaces
      if (!tempString.empty()) {
        std::stringstream sline(tempString);
        for (int i=0; i<nValGiven;i++) {
          double temp;
          sline >> temp;
          config->nxVec.push_back(temp);
        }
      }
    }
    if (iter->first == "NYVEC") {
      tempString = trim(iter->second);
      int nValGiven = std::count(tempString.begin(), tempString.end(), ',') + 1;
      std::replace(tempString.begin(), tempString.end(), ',', ' ' ); // replace the comas with spaces
      if (!tempString.empty()) {
        std::stringstream sline(tempString);
        for (int i=0; i<nValGiven;i++) {
          double temp;
          sline >> temp;
          config->nyVec.push_back(temp);
        }
      }
    }
    if (iter->first == "NZFILTVEC") {
      tempString = trim(iter->second);
      int nValGiven = std::count(tempString.begin(), tempString.end(), ',') + 1;
      std::replace(tempString.begin(), tempString.end(), ',', ' ' ); // replace the comas with spaces
      if (!tempString.empty()) {
        std::stringstream sline(tempString);
        for (int i=0; i<nValGiven;i++) {
          double temp;
          sline >> temp;
          config->nzFiltVec.push_back(temp);
        }
      }
    }
  }  
  if (std::count(tempNC.begin(), tempNC.end(), ',') < config->nbt-1) { 
      std::cout << std::count(tempNC.begin(), tempNC.end(), ',')+1 << " values of NC given while " << config->nbt;
      std::cout << " chains will run! " << std::endl;
      std::cout << "Terminating..." << std::endl;
      exit(0);
  }
  int numberOfParameters = 0;
  if(config->swaves)
    numberOfParameters = 2*config->npu;
  else
    numberOfParameters = 2*config->npu;
  tempNC = trim(tempNC);
  std::replace(tempNC.begin(), tempNC.end(), ',', ' ' ); // replace the comas with spaces
  if (!tempNC.empty()) {
    std::stringstream sline(tempNC);
    for (int i=0; i<config->nbt;i++) {
      double temp;
      sline >> temp;
      if(temp > numberOfParameters) {
        std::cout << "Number of component specified (" << temp << ") for NC is too high (this is a run with "<< numberOfParameters;
        std::cout << " parameters)" << std::endl;
        std::cout << "Terminating..." << std::endl;
        exit(0);
      }
      config->nc.push_back(temp);
    }
  }
  if(config->mpiConfig.rank == 0)
    checkConfig(config);  // Check the config
  if (config->verbose1 && config->mpiConfig.rank == 0)
    std::cout << "  Configuration file read" << std::endl;
}

void build_temperature_ladders(Configuration* config)
// Build the temperature ladders (T[0]=1)
{
  double incr=0.;
  config->T.push_back(TMIN);
  incr=pow(config->tmax/config->T[0],1.0/(config->nbt-1));
  for (int i=1;i<config->nbt;i++)
    config->T.push_back(config->T[i-1]*incr);
}

void loadArrivalTimes(Configuration* config)
// Load the first arrival times for P and S waves contained in the file tObs4096.txt
{
  std::string pick(config->filesDir+config->name_of_times_file);
  std::ifstream file(pick.c_str()); // First arrival times P and S
  if(file == NULL) {
    std::cout << "IMPOSSIBLE TO OPEN "+config->filesDir+config->name_of_times_file << std::endl;
    exit(0);
  }
  std::string line;
  while(getline(file, line)) {
    line = trim(line);
    if (!line.empty()) {
      std::stringstream sline(line);
      double tp,ts;
      sline >> tp >> ts;
      config->data.times.timesP.push_back(tp);
      config->data.times.timesS.push_back(ts);
    }
  }
  file.close();
  if (config->verbose1 && config->mpiConfig.rank == 0)
    std::cout << "  Arrival times loaded" << std::endl;
}

void loadStationsPositions(Configuration* config)
// Load the receivers positions contained in the file coordStats.txt
{
  std::string coordStations(config->filesDir+config->name_of_stations_file);
  std::ifstream file(coordStations.c_str()); // First arrival times P and S
  if(file == NULL) {
    std::cout << "IMPOSSIBLE TO OPEN "+config->filesDir+config->name_of_stations_file << std::endl;
    exit(0);
  }
  Coordinate X;
  std::string line;
  while(getline(file, line)) {
    line = trim(line);
    if (!line.empty()) {
      std::stringstream sline(line);
      sline >> X.x >> X.y >> X.z ;
      config->data.coordStations.push_back(X);
    }
  }
  file.close();
  if (config->verbose1 && config->mpiConfig.rank == 0)
    std::cout << "  Receivers positions loaded" << std::endl;
}

void loadShotsPositions(Configuration* config)
// Load the sources positions contained in the file coordShots.txt
{
  std::string coordShots(config->filesDir+config->name_of_shots_file);
  std::ifstream file(coordShots.c_str()); // First arrival times P and S
  if(file == NULL) {
    std::cout << "IMPOSSIBLE TO OPEN "+config->filesDir+config->name_of_shots_file << std::endl;
    exit(0);
  }
  Coordinate X;
  std::string line;
  while(getline(file, line)) {
    line = trim(line);
    if (!line.empty()) {
      std::stringstream sline(line);
      sline >> X.x >> X.y >> X.z ;
      config->data.coordShots.push_back(X);
    }
  }
  file.close();
  
  #ifdef PAR 
    int maxProcNeeded = 0;
    if(config->swaves)
      maxProcNeeded =  2*(int)config->data.coordShots.size();
    else
      maxProcNeeded = (int)config->data.coordShots.size();  
    if (config->mpiConfig.nb_process > maxProcNeeded) {
      if(config->mpiConfig.rank == 0) {
        std::cout << "You are using too much CPUs, " << maxProcNeeded << " is enough !" << std::endl;
        std::cout << "Terminating.." << std::endl;
      }
      exit(0);
    }
  #endif
  
  if (config->verbose1 && config->mpiConfig.rank == 0)
    std::cout << "  Shots positions loaded" << std::endl;
}

void loadPriorFeatures(Configuration* config)
// Load prior features (min, max, amplitude of variation) contained in the file priorFeatures.txt. Calculate prior energy.
// TODO write that if useAllWavelets
{
  std::string priorFeatures(config->filesDir+config->name_of_prior_features_file);
  std::ifstream file(priorFeatures.c_str()); // First arrival times P and S
  if(file == NULL) {
    std::cout << "IMPOSSIBLE TO OPEN "+config->filesDir+config->name_of_prior_features_file << std::endl;
    exit(0);
  }
  std::string line;
  while(getline(file, line)) {
    line = trim(line);
    if (!line.empty()) {
      std::stringstream sline(line);
      int idx;  // The params[i] refer to the wavelet coefficient number idx of the profile
      double mini,maxi;
      if (config->waveletParameterization)
        sline >> idx >> mini >> maxi ;
      else
        sline >> mini >> maxi ;
      config->data.minParameters[0].push_back(mini);
      config->data.maxParameters[0].push_back(maxi);
      if (config->waveletParameterization)
        config->data.indexParameters[0].push_back(idx);
    }
  }
  file.close();

  double Ep=1.;
  for(int i=0;i<(int)config->data.minParameters[0].size();i++) { // Loop on the parameters, to calculate the prior probability and computed the non varying parameters
    if(fabs(config->data.maxParameters[0][i]-config->data.minParameters[0][i]) > TINYVAL)
      Ep*=fabs(config->data.maxParameters[0][i]-config->data.minParameters[0][i]);
    else
      config->data.staticParameters.push_back(i); // This parameters are static (minParam=maxParam : they are actually not parameters...)
  }
  config->data.Ep=log(Ep);
  if (config->verbose1 && config->mpiConfig.rank == 0)
    std::cout << "  Prior features loaded" << std::endl;
}

void findOptimumFirstGuessParameterization(Configuration* config)
// From the first guess velocity profile(s) and given the parameterization scheme (store in config), determine the best set of
// coefficients describing it. Store it in config->coeffsP[wavelet] (and in config->coeffsS[wavelet]). Then store the corresponding parameterized profile
// in config->data.filtFirstGuessP[wavelet] (and in config->data.filtFirstGuessP[wavelet]).
{
  if (config->verbose1 && config->mpiConfig.rank == 0)
    std::cout << "Finding optimum first guess parameterization..." << std::endl;
  if (config->waveletParameterization) {
    for (int wavelet = 0; wavelet < config->nWavelets; wavelet++) { // Loop on the wavelets used
      // ******** DISCRETE WAVELET TRANSFORM IMPLEMENTATION*********"
      // Performing Non Linear Approximation by using only config->npu most significant coefficients
      // Coefficients in config->coeffsP[wavelet] (or config->coeffsS[wavelet]) are stored as following
      // config->coeffsP[wavelet] =[Appx(J-1) Detail(J-1) Detail(J-2) .... Detail(0)]
      dwt(config->data.firstGuessP, config->ndwts, config->listOfWavelets[wavelet], config->coeffsP[wavelet],config->flagP[wavelet], config->lengthP[wavelet]); // Performs J-Level DWT
//      if (config->verbose2 && config->mpiConfig.rank == 0) {
//        std::cout << "Wavelet: " << config->listOfWavelets[wavelet] << std::endl;
//        std::cout << "CoeffsP[" << config->listOfWavelets[wavelet] << "] : ";
//        for (std::vector<double>::const_iterator j = config->coeffsP[wavelet].begin(); j != config->coeffsP[wavelet].end(); ++j)
//          std::cout << *j << ' ';
//        std::cout << std::endl;
//      }
      if (config->keep_first_values != 1)
        keepNsignificantValues(&config->coeffsP[wavelet],config->npu); // Keep the N biggest absolute values in vector. Put the others to 0
      else
        keepNfirstValues(&config->coeffsP[wavelet],config->npu); // Keep the N biggest absolute values in vector. Put the others to 0
//      if (config->verbose2 && config->mpiConfig.rank == 0) {
//        std::cout << "CoeffsP[" << config->listOfWavelets[wavelet] << "] : ";
//        for (std::vector<double>::const_iterator j = config->coeffsP[wavelet].begin(); j != config->coeffsP[wavelet].end(); ++j)
//          std::cout << *j << ' ';
//        std::cout << std::endl << std::endl;
//      }
      idwt(config->coeffsP[wavelet],config->flagP[wavelet],config->listOfWavelets[wavelet],config->data.filtFirstGuessP[wavelet],config->lengthP[wavelet]);  // Performs IDWT with approximated coefficients
      if (config->swaves) {
        // ******** DISCRETE WAVELET TRANSFORM IMPLEMENTATION*********"
        dwt(config->data.firstGuessS, config->ndwts, config->listOfWavelets[wavelet], config->coeffsS[wavelet],config->flagS[wavelet], config->lengthS[wavelet]); // Performs J-Level DWT
//        if (config->verbose2 && config->mpiConfig.rank == 0) {
//          std::cout << "CoeffsS[" << config->listOfWavelets[wavelet] << "] : ";
//          for (std::vector<double>::const_iterator j = config->coeffsS[wavelet].begin(); j != config->coeffsS[wavelet].end(); ++j)
//            std::cout << *j << ' ';
//          std::cout << std::endl;
//        }
        if (config->keep_first_values != 1)
          keepNsignificantValues(&config->coeffsS[wavelet],config->npu); // Keep the N biggest absolute values in vector. Put the others to 0
        else
          keepNfirstValues(&config->coeffsS[wavelet],config->npu); // Keep the N biggest absolute values in vector. Put the others to 0
//        if (config->verbose2 && config->mpiConfig.rank == 0) {
//          std::cout << "CoeffsS[" << config->listOfWavelets[wavelet] << "] : ";
//          for (std::vector<double>::const_iterator j = config->coeffsS[wavelet].begin(); j != config->coeffsS[wavelet].end(); ++j)
//            std::cout << *j << ' ';
//          std::cout << std::endl << std::endl;
//        }
        idwt(config->coeffsS[wavelet],config->flagS[wavelet],config->listOfWavelets[wavelet],config->data.filtFirstGuessS[wavelet],config->lengthS[wavelet]);  // Performs IDWT with approximated coefficients
      }
    }
  }
  else { // Layer based parameterization
                                                                // zbottom z0   z1  ...    zNPU-2  ztop
    double dzNPU = (config->zbottom - config->ztop)/config->npu;//    |    o    o   ...       o    |     --> zbottom - ztop = (NPU-2)dz+2dz
                                                                //     <dz> <dz>               <dz>
    std::vector<double> zLayers,optimumZlayers,vp,vs;
    // Compute the first point :
    optimumZlayers.push_back(config->ztop+dzNPU);
    config->indices.push_back(closest(config->data.zp,optimumZlayers.back()));
    zLayers.push_back(config->data.zp[config->indices[0]]); // Value closest of z0 in zp
    // Construct vTemp from subrange in firstGuessP :
    std::vector<double> vTemp(config->data.firstGuessP.begin(), config->data.firstGuessP.begin() + config->indices[0]); // Extract the elements corresponding to first layer
    double meanValue = std::accumulate(vTemp.begin(), vTemp.end(), 0.0)/vTemp.size(); // Compute their mean value
    config->coeffsP[0].push_back(meanValue); // Store it in coeffsP[0]
    // Build z values for each layer
    for(int i=1;i<config->npu-1;i++) {
      optimumZlayers.push_back(optimumZlayers[i-1]+dzNPU);
      config->indices.push_back(closest(config->data.zp,optimumZlayers.back()));
      zLayers.push_back(config->data.zp[config->indices.back()]);
      // Construct vTemp2 from subrange in firstGuessP :
      std::vector<double> vTemp2(config->data.firstGuessP.begin()+config->indices[i-1], config->data.firstGuessP.begin() + config->indices[i] + 1); // Extract the elements corresponding to ith layer
      meanValue = std::accumulate(vTemp2.begin(), vTemp2.end(), 0.0)/vTemp2.size(); // Compute their mean value
      config->coeffsP[0].push_back(meanValue); // Store it in coeffsP[0]
    }
    // Construct vTemp from subrange in firstGuessP :
    std::vector<double> vTemp3(config->data.firstGuessP.begin()+config->indices[config->npu-2], config->data.firstGuessP.end()); // Extract the elements corresponding to ith layer
    meanValue = std::accumulate(vTemp3.begin(), vTemp3.end(), 0.0)/vTemp3.size(); // Compute their mean value
    config->coeffsP[0].push_back(meanValue); // Store it in coeffsP[0]
    for(unsigned int i=0;i<zLayers.size();i++)
      config->zStepCenter.push_back(zLayers[i] - dzNPU/2.0);
    config->zStepCenter.push_back(zLayers.back()+dzNPU/2.0);
    // Store the filtered profile:
    InverseLayerTransform(&config->data.filtFirstGuessP[0],&(config->coeffsP[0]), config, false);
    if (config->verbose1 && config->mpiConfig.rank == 0) {
      for(unsigned int i=0;i<config->coeffsP[0].size();i++)
        std::cout << "zStepCenter[" << i << "]: " << config->zStepCenter[i] << " config->coeffsP[0][" << i << "] :" << config->coeffsP[0][i] << std::endl;
    }
    if (config->swaves) {
      // Construct vTemp from subrange in firstGuessS :
      std::vector<double> vTemp4(config->data.firstGuessS.begin(), config->data.firstGuessS.begin() + config->indices[0]); // Extract the elements corresponding to first layer
      double meanValue = std::accumulate(vTemp4.begin(), vTemp4.end(), 0.0)/vTemp4.size(); // Compute their mean value
      config->coeffsS[0].push_back(meanValue); // Store it in coeffsS[0]
      // Build z values for each layer
      for(int i=1;i<config->npu-1;i++) {
        // Construct vTemp from subrange in firstGuessS :
        std::vector<double> vTemp5(config->data.firstGuessS.begin()+config->indices[i-1], config->data.firstGuessS.begin() + config->indices[i] + 1); // Extract the elements corresponding to ith layer
        meanValue = std::accumulate(vTemp5.begin(), vTemp5.end(), 0.0)/vTemp5.size(); // Compute their mean value
        config->coeffsS[0].push_back(meanValue); // Store it in coeffsS[0]
      }
      // Construct vTemp from subrange in firstGuessS :
      std::vector<double> vTemp6(config->data.firstGuessS.begin()+config->indices[config->npu-2], config->data.firstGuessS.end()); // Extract the elements corresponding to ith layer
      meanValue = std::accumulate(vTemp6.begin(), vTemp6.end(), 0.0)/vTemp6.size(); // Compute their mean value
      config->coeffsS[0].push_back(meanValue); // Store it in coeffsS[0]
      if (config->verbose1 && config->mpiConfig.rank == 0) {
        for(unsigned int i=0;i<config->coeffsS[0].size();i++)
          std::cout << "zStepCenter[" << i << "]: " << config->zStepCenter[i] << " config->coeffsS[0][" << i << "] :" << config->coeffsS[0][i] << std::endl;
      }
      // Store the filtered profile:
      InverseLayerTransform(&config->data.filtFirstGuessS[0],&(config->coeffsS[0]), config, true);
    }
  }
  if (config->verbose1 && config->mpiConfig.rank == 0)
    std::cout << "Done!" << std::endl;
}

void loadFirstGuesses(Configuration* config)
// Load first guesses in the files falseFirstGuessP4096.txt, falseFirstGuessS4096.txt. Filter them and store the results in config. 
// Must be done after loading priorFeatures.txt 
{
  std::string firstGuessP(config->filesDir+config->name_of_first_guess_P_file);
  std::ifstream fileP(firstGuessP.c_str());
  if(fileP == NULL) {
    std::cout << "IMPOSSIBLE TO OPEN "+config->filesDir+config->name_of_first_guess_P_file << std::endl;
    exit(0);
  }
  std::string lineP;
  while(getline(fileP, lineP)) {
    lineP = trim(lineP);
    if (!lineP.empty()) {
      std::stringstream slineP(lineP);
      double zpFirstGuess,value;
      slineP >> zpFirstGuess >> value ;
      config->data.firstGuessP.push_back(value);
      config->data.zpFirstGuess.push_back(zpFirstGuess);
    }
  }
  config->data.zp = config->data.zpFirstGuess; // In fact these two vectors are equal even for an analytical run but we need to declare it.
  if (config->data.zp[1] > config->data.zp[0])
    config->data.dz = config->data.zp[1]-config->data.zp[0];
  else {
    std::cout << "By now z must point downwards. Terminating..." << std::endl;
    exit(0);
  } 
  if(config->swaves) {
    std::string firstGuessS(config->filesDir+config->name_of_first_guess_S_file);
    std::ifstream fileS(firstGuessS.c_str());
    if(fileS == NULL) {
      std::cout << "IMPOSSIBLE TO OPEN "+config->filesDir+config->name_of_first_guess_S_file << std::endl;
      exit(0);
    }
    std::string lineS;
    int i=0;
    while(getline(fileS, lineS)) {
      lineS = trim(lineS);
      if (!lineS.empty()) {
        std::stringstream slineS(lineS);
        double zpFirstGuess,value;
        slineS >> zpFirstGuess >> value ;
        if (fabs(zpFirstGuess-config->data.zpFirstGuess[i]) > TINYVAL) {
          std::cout << "Depths read in "+config->filesDir+config->name_of_first_guess_S_file+" are not consistent with those in : ";
          std::cout << config->filesDir+config->name_of_first_guess_P_file+"! Terminating..." << std::endl;
          exit(0);
        }
        config->data.firstGuessS.push_back(value);
      }
      i++;
    }
    fileS.close();
    if (config->data.firstGuessS.size() != config->data.firstGuessP.size()) {
      std::cout << "First guess P and S waves velocities don't have the same number of points. Terminating..." << std::endl;
      exit(0);
    }
  }
  // Build z : it will contains the z positions of the cells borders
  config->data.z.push_back(config->data.zp[0]-config->data.dz/2.0); // Compute the first point
  if (fabs(config->data.z[0]) < TINYVAL)
    config->data.z[0] = TINYVAL;
  for(int i=1;i<(int)config->data.zp.size()+1;i++)
    config->data.z.push_back(config->data.z[i-1]+config->data.dz);
  fileP.close();
  config->data.nz=config->data.z.size();
  config->ztop=config->data.zp[0];
  config->zbottom=config->data.zp.back();
  
  if (! config->calculateTimesForFirstGuess) {
    findOptimumFirstGuessParameterization(config); // From the first guess velocity profile(s) and given the parameterization scheme
    // (store in config), determine the best set of coefficients describing it. Store it in config->coeffsP[wavelet] (and config->coeffsS[wavelet]).
    // Then store the corresponding parameterized profile in config->data.filtFirstGuessP[wavelet] (config->data.filtFirstGuessP[wavelet]).
    if (config->mpiConfig.rank == 0) {    // (In case of parallel implementation just one process has to create files)
      for (int wavelet = 0; wavelet < config->nWavelets; wavelet++) { // Loop on the wavelets used
        write_two_columns_file(&config->data.zp,&config->data.filtFirstGuessP[wavelet], config->outputDir+"filteredFirstGuessP."+config->listOfWavelets[wavelet]+"."+config->code+".dat");
        if (config->swaves)
          write_two_columns_file(&config->data.zp,&config->data.filtFirstGuessS[wavelet], config->outputDir+"filteredFirstGuessS."+config->listOfWavelets[wavelet]+"."+config->code+".dat");
      }
      if (config->verbose1)
        std::cout << "  First guess file loaded (and optimum filtering performed)" << std::endl;
    }
  }
}

void loadRealProfiles(Configuration* config)
// Load the real profiles from NAME_OF_REAL_PROFILE_FILE_P (and NAME_OF_REAL_PROFILE_FILE_S)
{
  std::string realP(config->filesDir+config->name_of_real_profile_P);
  std::ifstream fileP(realP.c_str());
  std::vector<double> zp;
  if(fileP == NULL) {
    std::cout << "IMPOSSIBLE TO OPEN "+config->filesDir+config->name_of_real_profile_P << std::endl;
    exit(0);
  }
  std::string lineP;
  while(getline(fileP, lineP)) {
    lineP = trim(lineP);
    if (!lineP.empty()) {
      std::stringstream slineP(lineP);
      double zpReal,velP;
      slineP >> zpReal >> velP ;
      config->data.realP.push_back(velP);
      zp.push_back(zpReal);
    }
  }
  fileP.close();
  if(zp.size() == config->data.zp.size()) {
    for (int i=0;i<(int)config->data.zp.size();i++) {
      if (fabs(zp[i]-config->data.zp[i]) > TINYVAL) {
        std::cout << "The depth from P wave velocity real profile are not consistent with those from P wave velocity first guess.";
        std::cout << "Terminating..." << std::endl;
        exit(0);
      }
    }
  }
  else {
    std::cout << "The depth from P wave velocity real profile are not consistent with those from P wave velocity first guess.";
    std::cout << "Terminating..." << std::endl;
    exit(0);
  }
  
  if(config->swaves) {
    std::string realS(config->filesDir+config->name_of_real_profile_S);
    std::ifstream fileS(realS.c_str());
    if(fileS == NULL) {
      std::cout << "IMPOSSIBLE TO OPEN "+config->filesDir+config->name_of_real_profile_S << std::endl;
      exit(0);
    }
    std::string lineS;
    int i=0;
    while(getline(fileS, lineS)) {
      lineS = trim(lineS);
      if (!lineS.empty()) {
        std::stringstream slineS(lineS);
        double zpReal,velS;
        slineS >> zpReal >> velS ;
        if (fabs(zpReal-config->data.zp[i]) > TINYVAL) {
          std::cout << "Depths read in "+config->filesDir+config->name_of_real_profile_S+" are not consistent with those in : ";
          std::cout << config->filesDir+config->name_of_real_profile_P+"! Terminating..." << std::endl;
          exit(0);
        }
        config->data.realS.push_back(velS);
      }
      i++;
    }
    fileS.close();
    if (config->verbose1 && config->mpiConfig.rank == 0)
      std::cout << "  Real velocity profiles loaded" << std::endl;
  }
  else
    if (config->verbose1 && config->mpiConfig.rank == 0)
      std::cout << "  Real velocity profile loaded" << std::endl;
}

void loadData(Configuration* config)
// Load the data (arrival times, first guesses, sources and receivers positions) contained in the files
{
  if (config->verbose1 && config->mpiConfig.rank == 0)
    std::cout << "Loading the data..." << std::endl;
  loadStationsPositions(config);
  loadShotsPositions(config);
  if (!config->buildPrior) {
    if (config->useAllWavelets) {
      std::cout << "For now: build Prior option is necessary when inversing wavelet type also..." << std::endl; // TODO
      exit(1);
    }
    else
      loadPriorFeatures(config);
  }
  loadFirstGuesses(config);
  if (config->analyticalRun) // If we perform an analytical run...
    loadRealProfiles(config); // ... load the real profiles
  else {// otherwise...
    if (! config->calculateTimesForFirstGuess) 
      loadArrivalTimes(config); // we load the arrival times
  }
  copyDataFiles(config); // Copy the data files used during the simulation in OUTPUT_FILES
  if (config->verbose1 && config->mpiConfig.rank == 0)
    std::cout << "Loading done !" << std::endl << std::endl;
}

void buildPrior(Configuration* config)
// If wavelet parameterization :
// Record the most significant coefficients of the filtered first guess profiles in config->data.indexParameters[wavelet]
// define their maximum variation ranges (config->data.minParameters[wavelet], config->data.maxParameters[wavelet]) during the rest of the algorithm.
// Calculate prior energy.
{
  std::vector<double> allMinParameters,allMaxParameters; // Will contain min,max for all wavelets
  std::vector<int> allIndexParameters; // Will contain index for all wavelets
  if (config->verbose1 && config->mpiConfig.rank == 0)
    std::cout << "Building the a priori space..." << std::endl;
  double Ep=0.0;
  for (int wavelet = 0; wavelet < config->nWavelets; wavelet++) { // Loop on the wavelets used
    for (unsigned int i = 0; i < config->coeffsP[wavelet].size(); i++) { // Loop on the P coeffs of first guess
      if (fabs(config->coeffsP[wavelet][i]) > 0) {
        config->data.maxParameters[wavelet].push_back((1.0+sign(config->coeffsP[wavelet][i])*config->A)*config->coeffsP[wavelet][i]);
        if (not config->waveletParameterization and (1.0-sign(config->coeffsP[wavelet][i])*config->A)*config->coeffsP[wavelet][i] < 0)
          config->data.minParameters[wavelet].push_back(TINYVAL); // If layer parameterization the parameters are velocities and can't be negative
        else
          config->data.minParameters[wavelet].push_back((1.0-sign(config->coeffsP[wavelet][i])*config->A)*config->coeffsP[wavelet][i]);
        if (config->waveletParameterization)
          config->data.indexParameters[wavelet].push_back((int)i);
      }
    }
    if (config->swaves) {
      for (unsigned int i = 0; i < config->coeffsS[wavelet].size(); i++) { // Loop on the S coeffs
        if (fabs(config->coeffsS[wavelet][i]) > 0) {
          config->data.maxParameters[wavelet].push_back((1.0+sign(config->coeffsS[wavelet][i])*config->A)*config->coeffsS[wavelet][i]);
          if (not config->waveletParameterization and (1.0-sign(config->coeffsS[wavelet][i])*config->A)*config->coeffsS[wavelet][i] < 0)
            config->data.minParameters[wavelet].push_back(TINYVAL); // If layer parameterization the parameters are velocities and can't be negative
          else
            config->data.minParameters[wavelet].push_back((1.0-sign(config->coeffsS[wavelet][i])*config->A)*config->coeffsS[wavelet][i]);
          if (config->waveletParameterization)
            config->data.indexParameters[wavelet].push_back((int)i);
        }
      }
    }

    // Concatenate into allMinParameters,allMaxParameters and allIndexParameters
    allMinParameters.insert(allMinParameters.end(), config->data.minParameters[wavelet].begin(), config->data.minParameters[wavelet].end() );
    allMaxParameters.insert(allMaxParameters.end(), config->data.maxParameters[wavelet].begin(), config->data.maxParameters[wavelet].end() );
    allIndexParameters.insert(allIndexParameters.end(), config->data.indexParameters[wavelet].begin(), config->data.indexParameters[wavelet].end() );

    for(int i=0;i<(int)config->data.minParameters[wavelet].size();i++) // Loop on the parameters to calculate the prior probability
      Ep+=log(fabs(config->data.maxParameters[wavelet][i]-config->data.minParameters[wavelet][i])); // TODO: verify that for layer based parameterization
    // Important! We use these vectors in InverseWaveletTransform to store the coefficients (that is stupid by the way...) so we
    // put it here back to zero TODO :
    for(int i=0;i<(int)config->coeffsP[wavelet].size();i++) 
      config->coeffsP[wavelet][i]=0.0;
    if(config->swaves) {
      for(int i=0;i<(int)config->coeffsS[wavelet].size();i++) 
        config->coeffsS[wavelet][i]=0.0;
    }
  }
  config->data.Ep=Ep;
  if (config->mpiConfig.rank == 0) {
    if (config->waveletParameterization)
      write_three_columns_file(&allIndexParameters,&allMinParameters, &allMaxParameters,config->outputDir+"priorFeatures."+config->code+".dat");
    else
      write_two_columns_file(&allMinParameters, &allMaxParameters,config->outputDir+"priorFeatures."+config->code+".dat");
  }
  if (config->verbose1 && config->mpiConfig.rank == 0)
  std::cout << "Done ! (energy of the prior : " << config->data.Ep << ")" << std::endl << std::endl;
}

void generate_profiles_from_prior(Configuration* config)  
// Generate config.nPriorProfiles profiles from a priori space in config.outputDir/priorProfilesXXX/
{
  if (config->verbose1 && config->mpiConfig.rank == 0 && config->nPriorProfiles > 0) {
    std::cout << "Storing " << config->nPriorProfiles << " profiles from the prior space in ";
    std::cout <<  config->outputDir+"priorProfiles"+config->code+"/" << std::endl;
  }
  for(int i=0;i<config->nPriorProfiles;i++) {
    std::vector<double> params;
    std::vector<double> priorProfileP, priorProfileS;
    if (config->verbose2 && config->mpiConfig.rank == 0)
      std::cout << "  Profile " << i+1 << " :" << std::endl;
    int wavelet = rand() % config->nWavelets; // Draw wavelet in the range 0 to config->nWavelets-1
    if (config->verbose2 && config->mpiConfig.rank == 0 && config->useAllWavelets && config->waveletParameterization)
      std::cout<< "  Wavelet drawn: " << config->listOfWavelets[wavelet] << "..." << std::endl;
    for(int j=0;j<(int)config->data.minParameters[wavelet].size();j++) {  // Loop on the number of parameters that we will modify
      params.push_back(Uniform(config->data.minParameters[wavelet][j],config->data.maxParameters[wavelet][j]));
      if (config->verbose2 && config->mpiConfig.rank == 0) {
        std::cout<< "  params["<< j <<"] : " << params[j] << "   ";
        std::cout<< "(Uniform between : "<< config->data.minParameters[wavelet][j] <<" and " << config->data.maxParameters[wavelet][j] << ") " << std::endl;
      }
    }
    if (config->verbose2 && config->mpiConfig.rank == 0)
      std::cout << std::endl;
    if (config->waveletParameterization)
      InverseWaveletTransform(&priorProfileP,&params,config,wavelet,false); // Calculate the P waves velocity profile corresponding to the wavelet coefficients
    else
      InverseLayerTransform(&priorProfileP,&params,config,false); // Calculate the P waves velocity profile corresponding to the layer coefficients
    if(config->swaves) {
      if (config->waveletParameterization)
        InverseWaveletTransform(&priorProfileS,&params,config,wavelet,true); // Calculate the S waves velocity profile corresponding to the wavelet coefficients
      else
        InverseLayerTransform(&priorProfileS,&params,config,true); // Calculate the S waves velocity profile corresponding to the wavelet coefficients
    }
    std::string mkdir_command="mkdir -p "+config->outputDir+"priorProfiles"+config->code+"/";
    int status;
    if (config->mpiConfig.rank == 0) {
      status = system(mkdir_command.c_str()); // TODO : This is a little bit dirty but there is no simple equivalent to "mkdir -p"
      if(status != 0) {
        std::cout << "Error while creating "+config->outputDir+"priorProfiles"+config->code+"/"+" (status = "<< status << ")" << std::endl;
        exit(1);
      }
    }
    std::ostringstream ii;   // Store i as a string
    ii << i;
    if (config->mpiConfig.rank == 0) {
      write_two_columns_file(&config->data.zp, &priorProfileP, config->outputDir+"priorProfiles"+config->code+"/"+"priorProfileP."+config->code+"."+ii.str()+".dat");
      if(config->swaves)
        write_two_columns_file(&config->data.zp, &priorProfileS, config->outputDir+"priorProfiles"+config->code+"/"+"priorProfileS."+config->code+"."+ii.str()+".dat");
    }
  }
  if (config->verbose1 && config->mpiConfig.rank == 0 && config->nPriorProfiles > 0)
    std::cout << "Done ! " << std::endl << std::endl;
}

void designSmallGrid(Configuration* config)
// Calculate optimum nx,ny,nzFilt,dx,dy,dzFilt to reduce computation time.
{
  if(config->verbose2 && config->mpiConfig.rank == 0) 
    std::cout << "Designing small grid..." << std::endl;
  // Calculate min and max values of coordinates given in coordStats and coordShots files
  double xmin,ymin,zmin,xmax,ymax,zmax; // Will contain min and max values of coordinates given in coordStats and coordShots files
  std::vector<double> xs,ys,zs; // Will contain all x and y coordinates
  for (int istat=0;istat<(int)config->data.coordStations.size();istat++) {
    xs.push_back(config->data.coordStations[istat].x);
    ys.push_back(config->data.coordStations[istat].y);
    zs.push_back(config->data.coordStations[istat].z);
  }
  for (int ishot=0;ishot<(int)config->data.coordShots.size();ishot++) {
    xs.push_back(config->data.coordShots[ishot].x);
    ys.push_back(config->data.coordShots[ishot].y);
    zs.push_back(config->data.coordShots[ishot].z);
  }
  xmax = *std::max_element(xs.begin(), xs.end());
  ymax = *std::max_element(ys.begin(), ys.end());
  zmax = *std::max_element(zs.begin(), zs.end());
  xmin = *std::min_element(xs.begin(), xs.end());
  ymin = *std::min_element(ys.begin(), ys.end());
  zmin = *std::min_element(zs.begin(), zs.end());

  double zminFirstGuess = *std::min_element(config->data.z.begin(), config->data.z.end()); // min values of z coordinates given in first guess file
  double zmaxFirstGuess = *std::max_element(config->data.z.begin(), config->data.z.end()); // max values of z coordinates given in first guess file

  if ((xmin == xmax) || (ymin == ymax) || (zmin == zmax)) { // (**)
    std::cout << "This is a 2D (or 1D) problem... Using FTeik 3D is not efficient, we should use FTeik2D (or a simple t=d/v ";
    std::cout << "if the problem is 1D) but it is not implemented for now. Terminating..." << std::endl;
    exit(0);
  }  
  if ((zmin < zminFirstGuess) || (zmax > zmaxFirstGuess)) {
    std::cout << "z minimum or maximum values given in "+config->filesDir+config->name_of_stations_file+" and ";
    std::cout << config->filesDir+config->name_of_shots_file+" are inconsistent with those described by ";
    std::cout << config->filesDir+config->name_of_first_guess_P_file+". Terminating..." << std::endl;
    std::cout << "zmin : " << zmin << " zminFirstGuess : "<< zminFirstGuess << " zmax : " << zmax << " zmaxFirstGuess : "<< zmaxFirstGuess << std::endl;
    exit(0);
  }
  
  double xminGrid = xmin - (xmax-xmin)*config->coordTol;
  double xmaxGrid = xmax + (xmax-xmin)*config->coordTol;
  double yminGrid = ymin - (ymax-ymin)*config->coordTol;
  double ymaxGrid = ymax + (ymax-ymin)*config->coordTol;
  
  config->data.xminGridGlob = xminGrid;
  config->data.yminGridGlob = yminGrid;
  config->data.zminGridGlob = config->data.z[0];
  if(config->verbose2 && config->mpiConfig.rank == 0) {
    std::cout << "  Min x coordinate : " << xmin << " Max x coordinate "<< xmax << " Min y coordinate : " << ymin << " Max y coordinate : "<< ymax << std::endl;
    std::cout << "  Min x of grid : " << xminGrid << " Max x of grid "<< xmaxGrid << " Min y of grid : " << yminGrid << " Max y of grid : "<< ymaxGrid << std::endl;
    std::cout << "  Min z : " << zmin << " Min z of grid : " << config->data.zminGridGlob << std::endl;
    std::cout << "  Max z : " << config->data.z.back() << " Min z of grid : " << config->data.zminGridGlob << std::endl << std::endl;
  }
  // Find min relatives intervals between two stations and/or receivers (expressed as a fraction of the maximum distances xmax and ymax)
  // Concatenate xstat with xshots -> xcoords, and ystat with yshots -> ycoords
  std::vector<double> dx,dy; // Will contain the relatives intervals
  for (int ix1=0;ix1<(int)xs.size()-1;ix1++) {
    for (int ix2=ix1+1;ix2<(int)xs.size();ix2++) {
      //std::cout << "(ix1,ix2) : (" << ix1 << "," << ix2 << ")" << std::endl;
      if (fabs(xs[ix1]-xs[ix2]) > TINYVAL) // We just look at distinct coordinates
        dx.push_back(fabs(xs[ix1]-xs[ix2]));
      if (fabs(ys[ix1]-ys[ix2]) > TINYVAL) // We just look at distinct coordinates
        dy.push_back(fabs(ys[ix1]-ys[ix2]));
    }
  }
  double dxmin = -1.0,dymin=-1.0;
  if (dx.size()>0) // If there is distinct coordinates. Otherwise the problem is 2D or 1D! This must not happen! See (**) above
    dxmin = *std::min_element(dx.begin(), dx.end()); // Minimum x distances between two stations and/or receivers
  if (dy.size()>0) // If there is distinct coordinates. Otherwise the problem is 2D or 1D! This must not happen! See (**) above
    dymin = *std::min_element(dy.begin(), dy.end()); // Minimum y distances between two stations and/or receivers
  
  if (config->calculateTimesForFirstGuess) {
    std::cout << "  Min x of grid : " << xminGrid << " Max x of grid "<< xmaxGrid << " Min y of grid : " << yminGrid << " Max y of grid : "<< ymaxGrid << std::endl;
    config->data.nx=config->nxDefault;
    config->data.ny=config->nyDefault;
    config->data.dx = (xmaxGrid-xminGrid)/((float)(config->data.nx-1));
    config->data.dy = (ymaxGrid-yminGrid)/((float)(config->data.ny-1));
    // nz, and dz are already known from first guess curves (function loadFirstGuesses)
  }
  else
    findOptimumGrid(xmaxGrid,ymaxGrid,dxmin,dymin,config); // if config->findOptimumGrid = 1 :
  // Filter a first guess curve iteratively and run the eikonal each time to determine nx,ny and nzFilt to use 
  // if config->findOptimumGrid = 0 use : 
  // nx = config->nxDefault
  // ny = config->nyDefault
  // nzFilt = config->nzfiltDefault
  //
  // --> Then calculate dx,dy and dzFilt
}

void calculateTimesForFirstGuess(Configuration* config)
// Runs the eikonal for first guess(es) curves
{

  if(config->mpiConfig.rank == 0) {
    std::cout << "For this run we will only calculate the times for the curve(s) given as :" << std::endl;
    std::cout << config->filesDir+config->name_of_first_guess_P_file << std::endl;
    if(config->swaves)
      std::cout << config->filesDir+config->name_of_first_guess_S_file << std::endl;
    if (config->resample)
      std::cout << "They will be resampled to NZ = " << config->nzfiltDefault << std::endl; 
    std::cout << "The sources/receivers given in :" << std::endl; 
    std::cout << config->filesDir+config->name_of_shots_file << std::endl;
    std::cout << config->filesDir+config->name_of_stations_file << std::endl;
  }
  
  VelocityModel velModel;    // Will contain the velocity model
  ArrivalTimes arrivalTimes; // Will contain the P and S waves arrival times at receivers positions
  std::vector<double> velP1D,velS1D;
  if (config->resample) {
    config->data.zFilt=linspace(config->data.z[0],config->data.z.back(),config->nzfiltDefault);
    // Build zFiltp : it will contains the real positions of the velocities vpFilt calculated after that
    for(int i=1;i<(int)config->data.zFilt.size();i++) // zFiltp give z in the middle of intervals
      config->data.zFiltp.push_back(config->data.zFilt[i-1]+(config->data.zFilt[i]-config->data.zFilt[i-1])/2.0);
    velP1D = downSampleProfile(config->data.firstGuessP,config->data.zp,config->data.zFiltp);
    if(config->mpiConfig.rank == 0)
      write_two_columns_file(&config->data.zFiltp,&velP1D, config->outputDir+"downSampledFirstGuessP."+config->code+".dat");
    double dzFilt=config->data.dz*((double)config->data.nz-1)/((double)config->nzfiltDefault-1.0);
    velModel.nx = config->data.nx; velModel.ny = config->data.ny; velModel.nz = config->nzfiltDefault;  
    velModel.dx = config->data.dx; velModel.dy = config->data.dy; velModel.dz = dzFilt;//fabs(config->data.z[1]-config->data.z[0]);
  }
  else {
    velModel.nx = config->data.nx; velModel.ny = config->data.ny; velModel.nz = config->data.nz;  
    velModel.dx = config->data.dx; velModel.dy = config->data.dy; velModel.dz = config->data.dz;//fabs(config->data.z[1]-config->data.z[0]);
  }
  calculateShiftForEachShot(&velModel,config);
  isThisConfigOk(&velModel,config);
  //    std::cout << "floor((config->data.coordShots[ishot].x-config->data.xminGridGlob)/config->data.dx) " << floor((config->data.coordShots[ishot].x-velModel.xmin[ishot])/velModel.dx) << std::endl;
  //  std::cout << "floor((config->data.coordShots[ishot].y-config->data.yminGridGlob)/config->data.dy) " << floor((config->data.coordShots[ishot].y-velModel.ymin[ishot])/velModel.dy) << std::endl;
  //  std::cout << "floor((config->data.coordShots[ishot].z-config->data.zminGridGlob)/config->data.dz) " << floor((config->data.coordShots[ishot].z-velModel.zmin[ishot])/velModel.dz) << std::endl;
  //  std::cout << "epsilonX " << epsilonX << std::endl;
  //  std::cout << "epsilonY " << epsilonY << std::endl;
  //  std::cout << "epsilonZ " << epsilonZ << std::endl;
  //  velModel.xmin = config->data.xminGrid; velModel.ymin = config->data.yminGrid; velModel.zmin = config->data.zminGrid;
  velModel.velP= new tab3d<double>(velModel.nz-1,velModel.nx-1,velModel.ny-1,-1.0); // Create a (nz-1,nx-1,ny-1) mesh and initialize every cell at -1
  // Copy the file content into the velocity model (extend the 1D profile to obtain a 3D profile) :
  if (config->resample)
    meshing(&velP1D,&velModel,false); // Extend this profile on the whole mesh
  else
    meshing(&config->data.firstGuessP,&velModel,false); // Extend this profile on the whole mesh
  
  std::cout << "Uses NX = " << velModel.nx << " NY = " << velModel.ny << " NZ = " << velModel.nz << std::endl;
  std::cout << "DX = " << velModel.dx << " DY = " << velModel.dy << " DZ = " << velModel.dz << std::endl;
  for(int shotNumber=0;shotNumber<(int)config->data.coordShots.size();shotNumber++) {
    std::cout << "Shot number " << shotNumber << std::endl;
    std::cout << "  xmin = " << velModel.xmin[shotNumber] << " ymin = " << velModel.ymin[shotNumber]  << " zmin = " << velModel.zmin[shotNumber]  << std::endl;
  }
  std::cout << "nSweeps = " << config->nSweeps << " epsin = " << config->epsin << std::endl << std::endl;
  
  tab3d<double> tt3dP(velModel.nz,velModel.nx,velModel.ny,1.0); // Will contain the P waves arrival times. Initialize every cell at 1.0 
  if(config->mpiConfig.rank == 0)
    std::cout << "Eikonal calculations for P waves... " << std::endl;
  for(int shotNumber=0;shotNumber<(int)config->data.coordShots.size();shotNumber++) {
    if(config->mpiConfig.rank == 0)
      std::cout << "  Shot Number " << shotNumber+1 << " on " << config->data.coordShots.size() << std::endl;  
    // Calculate the P waves travel times everywhere on the mesh (put them on tt3dP) for the shot number shotNumber :
    eikonal3d(&tt3dP,&velModel,config,shotNumber,false);
    for(int j=0;j<(int)config->data.coordStations.size();j++)
      arrivalTimes.timesP.push_back(getTime(&tt3dP,config->data.coordStations[j],&velModel,shotNumber));
     // arrivalTimes.timesP.push_back(getTime(&tt3dP,config->data.coordStations[j],&velModel,shotNumber)+0.01*shotNumber); // To create a synthetic case with false t0
  }
  delete velModel.velP;

  if(config->swaves) {
    velModel.velS= new tab3d<double>(velModel.nz-1,velModel.nx-1,velModel.ny-1,-1); // Create a (nz-1,nx-1,ny-1) mesh and initialize every cell at -1
    if (config->resample) {
      velS1D = downSampleProfile(config->data.firstGuessS,config->data.zp,config->data.zFiltp);
      if(config->mpiConfig.rank == 0)
        write_two_columns_file(&config->data.zFiltp,&velS1D, config->outputDir+"downSampledFirstGuessS."+config->code+".dat");
      meshing(&velS1D,&velModel,false); // Extend this profile on the whole mesh
    }
    else
      meshing(&config->data.firstGuessS,&velModel,true); // Extend this profile on the whole mesh
    tab3d<double> tt3dS(velModel.nz,velModel.nx,velModel.ny,1.0); // Will contain the S waves arrival times. Initialize every cell at 1.0 
    if(config->mpiConfig.rank == 0)
      std::cout << std::endl << "Eikonal calculations for S waves... " << std::endl;
    for(int shotNumber=0;shotNumber<(int)config->data.coordShots.size();shotNumber++) {
      if(config->mpiConfig.rank == 0)
        std::cout << "  Shot Number " << shotNumber+1 << " on " << config->data.coordShots.size() << std::endl;  
      // Calculate the S waves travel times everywhere on the mesh (put them on tt3dS) for the shot number shotNumber :
      eikonal3d(&tt3dS,&velModel,config,shotNumber,true);
      for(int j=0;j<(int)config->data.coordStations.size();j++)
        arrivalTimes.timesS.push_back(getTime(&tt3dS,config->data.coordStations[j],&velModel,shotNumber));
        //arrivalTimes.timesS.push_back(getTime(&tt3dS,config->data.coordStations[j],&velModel,shotNumber)+0.01*shotNumber); // To create a synthetic case with false t0
    }
    delete velModel.velS;
    if(config->mpiConfig.rank == 0)
      write_two_columns_file(&arrivalTimes.timesP,&arrivalTimes.timesS, config->outputDir+"calculatedTimes."+config->code+".dat");
  }
  else {
    if(config->mpiConfig.rank == 0)
      write_one_column_file(&arrivalTimes.timesP, config->outputDir+"calculatedTimes."+config->code+".dat");
  }
  if(config->mpiConfig.rank == 0) {
    std::cout << std::endl << "Done !" << std::endl;
    std::cout << "The calculated arrival times for first guess profiles has been written in : "+config->outputDir+"calculatedTimes."+config->code+".dat";
    std::cout << std::endl << std::endl;
  }
}

void findOptimumGrid(double xmaxGrid, double ymaxGrid, double dxmin, double dymin, Configuration* config)
// Filter a first guess curve iteratively and run the eikonal each time to determine nx,ny and nzFilt to use. 
// Then calculate dx,dy and dzFilt
// We just do that for P waves for now TODO : same for S waves as well
{
  int wavelet = 0;
  //int noXptdtsd = 0,noYptdtsd = 0; // Number Of Points To Describe The Smallest Distances
  int nx = 0, ny = 0, nzFilt = 0;
  float dx = 0.0, dy = 0.0, dzFilt = 0.0;
  float dxRef = 0.0, dyRef = 0.0, dzRef = 0.0;
  int shotNumber=config->shotNumberRef;
  std::vector<double> zFilt, zFiltp;
  std::vector<double> nxs,nys,nzs,dxs,dys,dzs,averageDiffs,percentOk; //noXptdtsdRecord,noYptdtsdRecord
  std::vector<int> oks;
   
  if(config->findOptimumGrid) {
    if(config->verbose1 && config->mpiConfig.rank == 0) {
      std::cout << "  -> Looking for optimum grid " << std::endl;
      if (config->useAllWavelets and config->waveletParameterization)
        std::cout << "    Wavelet: " << config->listOfWavelets[wavelet] << " will be used (first in the list)" << std::endl;
      std::cout << "    Some runs will be performed to find an optimum grid." << std::endl;
      std::cout << "    First, we construct references times with a precise grid. We extend the first guess P wave velocity profile on that grid." << std::endl;
      std::cout << "    (~"+formatBytes(sizeof(double)*config->nxref*config->nyref*config->data.nz+sizeof(double)*(config->nxref-1)*(config->nyref-1)*(config->data.nz-1))+" will be used... are you sure to have enough memory?)" << std::endl;
    }
    tab3d<double> refTt3dP(config->data.nz,config->nxref,config->nyref,-1.0); // Will contain the P waves arrival times. Initialize every cell at -1.0 
    VelocityModel refVelModel;    // Will contain the reference velocity model
    refVelModel.nx = config->nxref; refVelModel.ny = config->nyref; refVelModel.nz = config->data.nz;
    dxRef = (xmaxGrid-config->data.xminGridGlob)/(config->nxref-1);
    dyRef = (ymaxGrid-config->data.yminGridGlob)/(config->nyref-1);
    dzRef = config->data.dz;
    refVelModel.dx = dxRef; refVelModel.dy = dyRef; refVelModel.dz = dzRef;
    calculateShiftForEachShot(&refVelModel,config);
    isThisConfigOk(&refVelModel,config);
    if ((! refVelModel.OK) && (config->mpiConfig.rank == 0)) {
      std::cout << " Reference velocity model is not OK! Change NXREF, NYREF or COORD_TOL" << std::endl;
      std::cout << " Try to choose nx ~ ny ~ nz... but this can be prohibitive if the first guess curve contains two many points!" << std::endl;
      //  MPI_Finalize();                         // ...
      //  exit(0);                                // ... and we terminate!
    }
    //refVelModel.xmin = config->data.xminGrid; refVelModel.ymin = config->data.yminGrid; refVelModel.zmin = config->data.zminGrid;
    refVelModel.velP= new tab3d<double>(refVelModel.nz-1,refVelModel.nx-1,refVelModel.ny-1,-1.0); // Create a (nz-1,nx-1,ny-1) mesh and initialize every cell at -1
    meshing(&config->data.filtFirstGuessP[wavelet],&refVelModel,false); // Extend this profile on the whole mesh (it is the filteredFirstGuessP.XXX.dat)
    if(config->verbose1 && config->mpiConfig.rank == 0) {
      std::cout << "      Reference Eikonal computing for nx = config->nxref = " << config->nxref << " ny = config->nyref = " << config->nyref << " nz = ";
      std::cout << config->data.nz << "     (shot number " << shotNumber << ")... " << std::endl;
      std::cout << "      This could take some time..." << std::endl;
    }
    // Calculate the P waves travel times everywhere on the mesh (put them on tt3dP) for the shot number shotNumber :
    eikonal3d(&refTt3dP,&refVelModel,config,shotNumber,false); 
    delete refVelModel.velP;
    if(config->verbose1 && config->mpiConfig.rank == 0) {
      std::cout << "      Done !" << std::endl << std::endl;
      std::cout << "    Now we will compare the times from this large run to iteratively decimated version of P wave velocity profile... " << std::endl << std::endl;
    }    
    std::vector<int> nxVec,nyVec; // Nx,nyValues to test
    std::vector<int> nzFiltVec; // nz for filtered models

    nxVec=config->nxVec; //.push_back(4);nxVec.push_back(5);//nxVec.push_back(4);nxVec.push_back(5);
    nyVec=config->nyVec; //.push_back(4);nyVec.push_back(5);//nyVec.push_back(4);nyVec.push_back(5);
    nzFiltVec=config->nzFiltVec; //.push_back(100);nzFiltVec.push_back(150);//nzFiltVec.push_back(200);nzFiltVec.push_back(300);
    VelocityModel velModel;    // Will contain the velocity model
  
    for (int inzFilt=0;inzFilt<(int)nzFiltVec.size();inzFilt++) { 
      for (int inx=0;inx<(int)nxVec.size();inx++) {
        for (int iny=0;iny<(int)nyVec.size();iny++) {
          nzFilt = nzFiltVec[inzFilt];
        //  noXptdtsd = nxVec[inx];
        //  noYptdtsd = nyVec[iny];
        //  dx=dxmin/((float)noXptdtsd);
        //  dy=dymin/((float)noYptdtsd);
        //  nx=ceil((xmaxGrid-config->data.xminGridGlob)/dx);
        //  ny=ceil((ymaxGrid-config->data.yminGridGlob)/dy);
          nx = nxVec[inx];
          ny = nyVec[iny];
          dx = (xmaxGrid-config->data.xminGridGlob)/((float)(nx-1));
          dy = (ymaxGrid-config->data.yminGridGlob)/((float)(ny-1));
          dzFilt=config->data.dz*((double)config->data.nz-1)/((double)nzFilt-1.0);
          /*
          z[1]   -        zFilt[1]     -            ^
                dz                  dzFilt          |
          z[2]   -                                  |
                          zFilt[2]     -            |
          z[3]   -   --->                           |
                                                    |
          z[4]   -        zFilt[3]     -            | zmax-zmin = (nz-1)*dz = (nzFilt-1)*dzFilt   =>   dzFilt = dz*(nz-1)/(nzFilt-1)
                                                    |
          z[5]   -                                  |
                                                    |
           ...                ...                   |
                                                    |
          z[nz] -        zFilt[nzFilt] -            v
         
         !! Warning indices start at 0 (to nz-1) !! 
          */
          
          nxs.push_back(nx);
          nys.push_back(ny);
          nzs.push_back(nzFilt);
          dxs.push_back(dx);
          dys.push_back(dy);
          dzs.push_back(dzFilt);
          //noXptdtsdRecord.push_back(noXptdtsd);
          //noYptdtsdRecord.push_back(noYptdtsd);
          if(config->verbose2 && config->mpiConfig.rank == 0) {
            //std::cout << "    Number Of x Points To Describe The Smallest Distances : " << noXptdtsd << std::endl; 
            //std::cout << "    Number Of y Points To Describe The Smallest Distances : " << noYptdtsd << std::endl;
            std::cout << "    nx = " << nx << " ny = " << ny << " nzFilt = " << nzFilt;
            std::cout << "    dx = " << dx << " dy = " << dy << " dzFilt = " << dzFilt << std::endl;
            std::cout << "    (~"+formatBytes(sizeof(double)*nx*ny*nzFilt+sizeof(double)*(nx-1)*(ny-1)*(nzFilt-1))+" will be used... are you sure to have enough memory?)" << std::endl;
          }
          zFilt=linspace(config->data.z[0],config->data.z.back(),nzFilt);
          // Build zFiltp : it will contains the real positions of the velocities vpFilt calculated after that
          zFiltp.clear(); // Clear the vector leaving its size to 0
          for(int i=1;i<(int)zFilt.size();i++) // zFiltp give z in the middle of intervals
            zFiltp.push_back(zFilt[i-1]+(zFilt[i]-zFilt[i-1])/2.0);
          velModel.nx = nx; velModel.ny = ny; velModel.nz = nzFilt;
          velModel.dx = dx; velModel.dy = dy; velModel.dz = dzFilt;
          velModel.xmin.clear(); velModel.ymin.clear(); velModel.zmin.clear();
          calculateShiftForEachShot(&velModel,config);
          isThisConfigOk(&velModel,config);
          oks.push_back(velModel.OK);
          percentOk.push_back(velModel.percentOk);
          //velModel.xmin = config->data.xminGrid; velModel.ymin = config->data.yminGrid; velModel.zmin = config->data.zminGrid;
          std::vector<double> velP1D = downSampleProfile(config->data.filtFirstGuessP[wavelet],config->data.zp,zFiltp);
          // velP1D Contains the down sampled velocity model !! warning velP will not give the velocity at zFilt but between the points. It will have a size nzFilt-1
          velModel.velP= new tab3d<double>(velModel.nz-1,velModel.nx-1,velModel.ny-1,-1.0); // Create a (nz-1,nx-1,ny-1) mesh and initialize every cell at -1 
            // Copy the file content into the velocity model (extend the 1D profile to obtain a 3D profile).
          meshing(&velP1D,&velModel,false); // Extend this profile on the whole mesh
          tab3d<double> tt3dP(nzFilt,nx,ny,-1.0); // Will contain the P waves arrival times. Initialize every cell at -1.0 
          if(config->verbose2 && config->mpiConfig.rank == 0)
            std::cout << "      Eikonal computing for this model... " << std::endl;
          // Calculate the P waves travel times everywhere on the mesh (put them on tt3dP) for the shot number shotNumber :
          eikonal3d(&tt3dP,&velModel,config,shotNumber,false); 
          // Here compare tt with a big run
          if(config->verbose2 && config->mpiConfig.rank == 0)
            std::cout << "      Done !" << std::endl;
          delete velModel.velP;
          double averageDiff = 0;
          for (int k=0; k<(int)config->data.coordStations.size();k++) {
           // std::cout << "  tt3dP[" << k << "] : "  << getTime(&tt3dP,config->data.coordStations[k],&velModel,shotNumber);
           // std::cout << "  refTt3dP[ " << k << "] : " << getTime(&refTt3dP,config->data.coordStations[k],&refVelModel,shotNumber);
           // std::cout << "  diff : " << fabs(getTime(&tt3dP,config->data.coordStations[k],&velModel,shotNumber) -getTime(&refTt3dP,config->data.coordStations[k],&refVelModel,shotNumber)) << std::endl;
            averageDiff += fabs(getTime(&tt3dP,config->data.coordStations[k],&velModel,shotNumber) -getTime(&refTt3dP,config->data.coordStations[k],&refVelModel,shotNumber));
          }
          averageDiff /= (double)config->data.coordStations.size();
          averageDiffs.push_back(averageDiff);
          if(config->verbose2 && config->mpiConfig.rank == 0)
            std::cout << "    Average differences in seconds at the receivers between big run and that run for shot number "<< shotNumber << " : "  << averageDiff << "  --> OK? " << velModel.OK << std::endl << std::endl;
        }
      }
    }
    
    std::vector<double> nxsOK,nysOK,nzsOK,dxsOK,dysOK,dzsOK,averageDiffsOK;
    for (int k=0; k<(int)oks.size();k++) {
      if (oks[k] == 1) {
        nxsOK.push_back(nxs[k]);
        nysOK.push_back(nys[k]);
        nzsOK.push_back(nzs[k]);
        dxsOK.push_back(dxs[k]);
        dysOK.push_back(dys[k]);
        dzsOK.push_back(dzs[k]);
        averageDiffsOK.push_back(averageDiffs[k]);
      }
    }
    
    if ((int)nxsOK.size() > 0) { // At least one of the configurations tested is ok
      int min_index = std::min_element(averageDiffsOK.begin(), averageDiffsOK.end()) - averageDiffsOK.begin(); // This gives the index of the minimum average value
      if(config->verbose1 && config->mpiConfig.rank == 0) {
        std::cout << "    The best grid configuration among those tested which are ok ("<< (int)averageDiffsOK.size() << " OK tested) seems to be :" << std::endl;
        std::cout << "      nx = " << nxsOK[min_index] << " ny = " << nysOK[min_index] << " nzFilt = " << nzsOK[min_index] << std::endl;
        std::cout << "      dx = " << dxsOK[min_index] << " dy = " << dysOK[min_index] << " dzFilt = " << dzsOK[min_index] << std::endl;
        std::cout << "      -> Average differences at the receivers between big run and that run for shot number "<< shotNumber << " : "  << averageDiffsOK[min_index] << std::endl << std::endl;
        std::cout << "    We keep this grid configuration !" << std::endl;
        std::cout << std::endl;
        std::cout << "    The other acceptable configuration are:" << std::endl;
        for (int k=0; k<(int)nxsOK.size();k++) {
          std::cout << "      nx = " << nxsOK[k] << " ny = " << nysOK[k] << " nzFilt = " << nzsOK[k] << std::endl;
          std::cout << "      dx = " << dxsOK[k] << " dy = " << dysOK[k] << " dzFilt = " << dzsOK[k] << std::endl;
          std::cout << "      Average Diff : " << averageDiffsOK[k] << std::endl;
          std::cout << std::endl;
        }
        std::cout << "Done!" << std::endl << std::endl;
      }
      config->data.nzFilt=nzs[min_index];
      config->data.nx=nxs[min_index];
      config->data.ny=nys[min_index];
      config->data.dzFilt=dzs[min_index];
      config->data.dx=dxs[min_index];
      config->data.dy=dys[min_index];
    }
    else {
      int min_index = std::min_element(averageDiffs.begin(), averageDiffs.end()) - averageDiffs.begin(); // This gives the index of the minimum average value
      int max_percent_index = std::max_element(percentOk.begin(), percentOk.end()) - percentOk.begin(); // This gives the index of the maximum percentOK
      if(config->verbose1 && config->mpiConfig.rank == 0) {
        std::cout << "    None configuration is OK!!" << std::endl;
        std::cout << "    The best grid configuration among those tested ("<< (int)averageDiffs.size() << " tested) seems to be :" << std::endl;
        std::cout << "      nx = " << nxs[max_percent_index] << " ny = " << nys[max_percent_index] << " nzFilt = " << nzs[max_percent_index] << std::endl;
        std::cout << "      dx = " << dxs[max_percent_index] << " dy = " << dys[max_percent_index] << " dzFilt = " << dzs[max_percent_index] << std::endl;
//        std::cout << "      nx = " << nxs[min_index] << " ny = " << nys[min_index] << " nzFilt = " << nzs[min_index] << std::endl;
//        std::cout << "      dx = " << dxs[min_index] << " dy = " << dys[min_index] << " dzFilt = " << dzs[min_index] << std::endl;
        std::cout << "      -> Average differences at the receivers between big run and that run for shot number "<< shotNumber << " : "  << averageDiffs[max_percent_index] << std::endl << std::endl;
        std::cout << "    We keep this grid configuration !" << std::endl;
        std::cout << "Done!" << std::endl << std::endl;
      }
      config->data.nzFilt=nzs[max_percent_index];
      config->data.nx=nxs[max_percent_index];
      config->data.ny=nys[max_percent_index];
      config->data.dzFilt=dzs[max_percent_index];
      config->data.dx=dxs[max_percent_index];
      config->data.dy=dys[max_percent_index];    
    }
  } // If config->findOptimumGrid == 0, we use values given in the configuration file :
  else {
    config->data.nzFilt=config->nzfiltDefault;
    config->data.dx=dxmin/((float)config->nxDefault);
    config->data.dy=dymin/((float)config->nyDefault);
    config->data.nx=ceil((xmaxGrid-config->data.xminGridGlob)/config->data.dx);
    config->data.ny=ceil((ymaxGrid-config->data.yminGridGlob)/config->data.dy);
    config->data.dzFilt=config->data.dz*((double)config->data.nz-1)/((double)config->data.nzFilt-1.0);
    if(config->verbose1 && config->mpiConfig.rank == 0) {
      std::cout << "  -> We use the grid configuration given :" << std::endl;
      //std::cout << "    Number Of x Points To Describe The Smallest Distances : " << config->nxDefault << std::endl; 
      //std::cout << "    Number Of y Points To Describe The Smallest Distances : " << config->nyDefault << std::endl;
      std::cout << "    nx = " << config->data.nx << " (= " << config->nxDefault << " ??)  ny = " << config->data.ny << " (= " << config->nyDefault << "??)  nzFilt = " << config->data.nzFilt << std::endl;
      std::cout << "    dx = " << config->data.dx<< " dy = " << config->data.dy << " dzFilt = " << config->data.dzFilt << std::endl;
    }
  }
  config->data.zFilt=linspace(config->data.z[0],config->data.z.back(),config->data.nzFilt);
  // Build zFiltp : it will contains the real positions of the velocities vpFilt calculated after that
  for(int i=1;i<(int)config->data.zFilt.size();i++) // zFiltp give z in the middle of intervals
    config->data.zFiltp.push_back(config->data.zFilt[i-1]+(config->data.zFilt[i]-config->data.zFilt[i-1])/2.0);
  std::vector<double> velP1D = downSampleProfile(config->data.filtFirstGuessP[wavelet],config->data.zp,config->data.zFiltp);
  if(config->mpiConfig.rank == 0)
    write_two_columns_file(&config->data.zFiltp,&velP1D, config->outputDir+"downSampledFiltredFirstGuessP."+config->code+".dat");
} 

void createDataset(Configuration* config)
// Create the arrival times from the real profile. TODO : parallelize this function
{
  if(config->swaves) { 
    if(config->verbose1 && config->mpiConfig.rank == 0) {
      std::cout << "This is an analytical run : we create the dataset from the real velocity profiles given... ";
      std::cout << "(files : " << config->filesDir+config->name_of_real_profile_P << " and " ;
      std::cout << config->filesDir+config->name_of_real_profile_S << ")"<< std::endl;
    }
  }
  else {
    if(config->verbose1 && config->mpiConfig.rank == 0) {
      std::cout << "This is an analytical run, we first begin by creating the dataset from the real velocity profile...";
      std::cout << "(file : " << config->filesDir+config->name_of_real_profile_P << ")" << std::endl << std::endl;
    }
  }
  VelocityModel velModel;    // Will contain the velocity model
  ArrivalTimes arrivalTimes; // Will contain the P and S waves arrival times at receivers positions
  if (config->data.z[1] > config->data.z[0])
    config->data.dz = config->data.z[1]-config->data.z[0]; // In this case (analytical) we don't use the default dz but we calculate it
  else {
    std::cout << "By now z must point downwards. Terminating..." << std::endl;
    exit(0);
  }
  config->data.nz = config->data.realP.size()+1; // We add one because there is on border point more than intervals (see README)
  velModel.nx = config->data.nx; velModel.ny = config->data.ny; velModel.nz = config->data.nz;  
  velModel.dx = config->data.dx; velModel.dy = config->data.dy; velModel.dz = config->data.dz;//fabs(config->data.z[1]-config->data.z[0]);
  calculateShiftForEachShot(&velModel,config);
  //velModel.xmin = config->data.xminGrid; velModel.ymin = config->data.yminGrid; velModel.zmin = config->data.zminGrid;
  velModel.velP= new tab3d<double>(velModel.nz-1,velModel.nx-1,velModel.ny-1,-1.0); // Create a (nz-1,nx-1,ny-1) mesh and initialize every cell at -1 
  // Copy the file content into the velocity model (extend the 1D profile to obtain a 3D profile).
  meshing(&config->data.realP,&velModel,false); // Extend this profile on the whole mesh
  tab3d<double> tt3dP(velModel.nz,velModel.nx,velModel.ny,1.0); // Will contain the P waves arrival times. Initialize every cell at 1.0 
  if(config->verbose1 && config->mpiConfig.rank == 0)
    std::cout << "Eikonal calculations for P waves... " << std::endl;
  for(int shotNumber=0;shotNumber<(int)config->data.coordShots.size();shotNumber++) 
  {
    if(config->verbose1 && config->mpiConfig.rank == 0)
      std::cout << "  Shot Number " << shotNumber+1 << " on " << config->data.coordShots.size() << std::endl;  
    // Calculate the P waves travel times everywhere on the mesh (put them on tt3dP) for the shot number shotNumber :
    eikonal3d(&tt3dP,&velModel,config,shotNumber,false); 
    for(int j=0;j<(int)config->data.coordStations.size();j++)
      arrivalTimes.timesP.push_back(getTime(&tt3dP,config->data.coordStations[j],&velModel,shotNumber));
      // arrivalTimes.timesP.push_back(getTime(&tt3dP,config->data.coordStations[j],&velModel,shotNumber)+0.01*shotNumber); // To create a synthetic case with false t0
  }
  delete velModel.velP;
  config->data.times.timesP = arrivalTimes.timesP;

  if(config->swaves) {
    velModel.velS= new tab3d<double>(velModel.nz-1,velModel.nx-1,velModel.ny-1,-1); // Create a (nz-1,nx-1,ny-1) mesh and initialize every cell at -1
    meshing(&config->data.realS,&velModel,true); // Extend this profile on the whole mesh
    tab3d<double> tt3dS(velModel.nz,velModel.nx,velModel.ny,1.0); // Will contain the S waves arrival times. Initialize every cell at 1.0 
    if(config->verbose1 && config->mpiConfig.rank == 0)
      std::cout << std::endl << "Eikonal calculations for S waves... " << std::endl;
    for(int shotNumber=0;shotNumber<(int)config->data.coordShots.size();shotNumber++) 
    {
      if(config->verbose1 && config->mpiConfig.rank == 0)
        std::cout << "  Shot Number " << shotNumber+1 << " on " << config->data.coordShots.size() << std::endl;  
       // Calculate the S waves travel times everywhere on the mesh (put them on tt3dS) for the shot number shotNumber :
      eikonal3d(&tt3dS,&velModel,config,shotNumber,true);
      for(int j=0;j<(int)config->data.coordStations.size();j++)
        arrivalTimes.timesS.push_back(getTime(&tt3dS,config->data.coordStations[j],&velModel,shotNumber));
        //arrivalTimes.timesS.push_back(getTime(&tt3dS,config->data.coordStations[j],&velModel,shotNumber)+0.01*shotNumber); // To create a synthetic case with false t0
    }
    delete velModel.velS;
    config->data.times.timesS = arrivalTimes.timesS;
    if(config->mpiConfig.rank == 0)
      write_two_columns_file(&config->data.times.timesP,&config->data.times.timesS, config->outputDir+"calculatedTimes."+config->code+".dat");
  }
  else {
    if(config->mpiConfig.rank == 0)
      write_one_column_file(&config->data.times.timesP, config->outputDir+"calculatedTimes."+config->code+".dat");
    if(config->verbose1 && config->mpiConfig.rank == 0) {
      std::cout << "Done !" << std::endl;
      std::cout << "The calculated arrival times has been written in : "+config->outputDir+"calculatedTimes."+config->code+".dat";
      std::cout << std::endl << std::endl;
    }
  }
}

State init_state(Configuration* config)
// Initialize a markov chain state
{
  State state;
  if (config->swaves)
    state.params=std::vector<double>(config->npu*2); // Create a vector that will contain the parameters describing the state
  else
    state.params=std::vector<double>(config->npu); // Create a vector that will contain the parameters describing the state
  int wavelet = rand() % config->nWavelets; // Draw wavelet in the range 0 to config->nWavelets -1
  if (config->verbose2 && config->mpiConfig.rank == 0 && config->useAllWavelets)
    std::cout<< "  Wavelet drawn: " << config->listOfWavelets[wavelet] << "..." << std::endl;
  state.wavelet=wavelet;
  for(int i=0;i<(int)state.params.size();i++) { // Loop on the parameters
    // Draw the parameter number i of the state in an uniform prior
    state.params[i]=Uniform(config->data.minParameters[state.wavelet][i],config->data.maxParameters[state.wavelet][i]);
    state.E=0;  // We need the value of the temperature to calculate the energy and so the chain in which it is,
    //for the moment we are just considering a single state out of a chain
    if (config->verbose2 && config->mpiConfig.rank == 0) {
      std::cout<< "Init state param["<< i <<"] : " << state.params[i] << "   ";
      std::cout<< "(Uniform between : "<< config->data.minParameters[state.wavelet][i] << " and " << config->data.maxParameters[state.wavelet][i] << ") " << std::endl;
    }
  }
  return state;
}

Run init_run(Configuration* config)
// Initialize the run, create each chain at each temperature and the first state of each one.
{
  Run run;
  std::vector<long double> line(config->nbt-1); // This will contain the first line of the matrix of importance weights
  for(int iz=0;iz<config->data.nzFilt-1;iz++) {
    run.averageP.push_back(0.0);
    run.varP.push_back(0.0);
    if (config->swaves) {
      run.averageS.push_back(0.0);
      run.varS.push_back(0.0);
      run.varVpVs.push_back(0.0);
    }
  }
  for(int i=0;i<config->nbt;i++) {              // Loop on the temperatures : i.e the number of chains
    if(config->verbose1 && config->mpiConfig.rank == 0)
      std::cout << "Initialization of chain " << i+1 << " on " << config->nbt << " (Temperature : "<< config->T[i] << ") " << std::endl;
    State initialState;                     // Create a new state
    initialState=init_state(config);        // And initialize it
    Chain* chain = new Chain();            // Create a pointer on a Chain
    chain->states.push_back(initialState); // Put the initialState on it
    chain->nc=config->nc[i];               // Copy from the config the nb of component to modify in a MH iteration for this chain
    chain->T=config->T[i];                 // Copy from the config the temperature of this chain
    chain->i=i;                            // Store the index of the chain
    // Initialize the statistics of the chain :
    chain->at=0;
    chain->rt=0;
    chain->od=0;
    chain->ps=0;
    chain->as=0;
    chain->rs=0;
    std::vector<double> emptyVectorOfDouble;
    for (int wavelet = 0; wavelet < config->nWavelets; wavelet++) // Loop on the wavelets used
      chain->deltaParameters.push_back(emptyVectorOfDouble);
    for (int wavelet = 0; wavelet < config->nWavelets; wavelet++) { // Loop on the wavelets used
      std::vector<double> L;
      for(int j=0;j<(int)chain->states[0].params.size();j++) {
        L.push_back(config->data.maxParameters[wavelet][j]-config->data.minParameters[wavelet][j]);
        // At T=Tmax deltaState[i]=L[i]/di
        // At T=1    deltaState[i]=L[i]/df
        double a = 0.;
        double b = 0.;
        if (fabs(config->tmax-1) > TINYVAL) { // false if prior iteration for example
          a = L[j]*(1./config->df-1./config->di)/(1.-config->tmax);
          b = L[j]*(1./config->di-config->tmax/config->df)/(1.-config->tmax);
        } 
        chain->deltaParameters[wavelet].push_back(a*chain->T+b); // The variation range of the parameters varies with the temperature
      }
    }
    chain->dotProductP = new mpfr_t[config->data.nzFilt-1];
    chain->dotProductVarP = new mpfr_t[config->data.nzFilt-1];
    chain->dotProductS = new mpfr_t[config->data.nzFilt-1];
    chain->dotProductVarS = new mpfr_t[config->data.nzFilt-1];
    chain->dotProductVarVpVs = new mpfr_t[config->data.nzFilt-1];

    run.chains.push_back(chain); // Add this chain to the run
    // We will keep all the profiles generated by the chains to calculate the quantiles :
    for(int iz=0;iz<config->data.nzFilt-1;iz++) {
      std::vector<double> values;
      run.chains[i]->profilesP.push_back(values); // The profiles will be filled by "energy" function (in functions.cpp)
      if (config->swaves)
        run.chains[i]->profilesS.push_back(values);
    }
    run.chains[i]->states[0].E=energy(&initialState,chain,config); // Compute the energy of the first state of the chain
    mpfr_set_default_prec(PREC); // Set default precision to PREC
    mpfr_init_set_d(run.chains[i]->sumOfweights,0.0,RND); // Set that also to 0;
    for(int iz=0;iz<config->data.nzFilt-1;iz++) {
      run.chains[i]->minP.push_back(run.chains[i]->profilesP[iz].back()); // When there is just one profile it is the minimum and the maximum :
      run.chains[i]->maxP.push_back(run.chains[i]->profilesP[iz].back());
      if (config->swaves) {
        run.chains[i]->minS.push_back(run.chains[i]->profilesS[iz].back());
        run.chains[i]->maxS.push_back(run.chains[i]->profilesS[iz].back());
      }
      run.chains[i]->averageP.push_back(run.chains[i]->profilesP[iz].back()); //At the beginning the average profile is the first one
      run.chains[i]->varP.push_back(0.);
      run.chains[i]->qInfP.push_back(0.);
      run.chains[i]->qSupP.push_back(0.);
      mpfr_init_set_d(run.chains[i]->dotProductP[iz],0.0,RND); // Set that also to 0;
      mpfr_init_set_d(run.chains[i]->dotProductVarP[iz],0.0,RND); // Set that also to 0;
      run.chains[i]->weightedAverageP.push_back(0.);
      run.chains[i]->weightedVarP.push_back(0.);
      if (config->swaves) {
        run.chains[i]->averageS.push_back(run.chains[i]->profilesS[iz].back());
        run.chains[i]->varS.push_back(0.);
        run.chains[i]->qInfS.push_back(0.);
        run.chains[i]->qSupS.push_back(0.);
        mpfr_init_set_d(run.chains[i]->dotProductS[iz],0.0,RND); // Set that also to 0;
        mpfr_init_set_d(run.chains[i]->dotProductVarS[iz],0.0,RND); // Set that also to 0;
        mpfr_init_set_d(run.chains[i]->dotProductVarVpVs[iz],0.0,RND); // Set that also to 0;
        run.chains[i]->weightedAverageS.push_back(0.);
        run.chains[i]->weightedVarS.push_back(0.);
        run.chains[i]->weightedVarVpVs.push_back(0.);
      }
    }

    if (i>0) {
      // Initialization of importance weights matrix and normalization coefficients.
      line[i-1]=(long double)expl(run.chains[i]->states[0].E *(1.-run.chains[i]->T/run.chains[i-1]->T)); // TODO change that (mpfr)
    }
    //if (i>0)
    //  line[i-1]=exp(-(run.chains[i-1]->states[0].E-run.chains[i]->states[0].E));
    //  line[i-1]=exp(run.chains[i-1]->states[0].E *(run.chains[i-1]->T/run.chains[i]->T -1));
    // Initialization of importance weights matrix and normalization coefficients.
    // SCI[0][l]=exp(Ep*T[i]*(1/T[i+1]-1/T[i]));
    // delete(chain);
  }
  run.minP=run.chains[0]->minP; // First we store the profile for chain 0 (we will compare it with the others)
  run.maxP=run.chains[0]->maxP;
  if (config->swaves) {
    run.minS=run.chains[0]->minS; // First we store the profile for chain 0 (we will compare it with the others)
    run.maxS=run.chains[0]->maxS;
  }
  for(int i=1;i<config->nbt;i++) { // Loop on the other chains
    for(int iz=0;iz<config->data.nzFilt-1;iz++) { // Loop on the depths
      if (run.minP[iz] > run.chains[i]->minP[iz])
        run.minP[iz] = run.chains[i]->minP[iz];
      if (run.maxP[iz] < run.chains[i]->maxP[iz])
        run.maxP[iz] = run.chains[i]->maxP[iz];
      if (config->swaves) {
        if (run.minS[iz] > run.chains[i]->minS[iz])
          run.minS[iz] = run.chains[i]->minS[iz];
        if (run.maxS[iz] < run.chains[i]->maxS[iz])
          run.maxS[iz] = run.chains[i]->maxS[iz];
      }
    }
  }
  run.bestE.push_back(run.chains[0]->states[0].E); // We do that just to put something in that vector (we will fill it with writeBestProfiles)
  run.idxE.push_back(0);
  run.chainBestE.push_back(0);

  run.SCI.push_back(line);
  return run;
}

