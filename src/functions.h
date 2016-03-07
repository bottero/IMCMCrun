#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "structures.h"

double getTime(tab3d<double>* tt3d,Coordinate coord, VelocityModel* velModel, int ishot);
// Give the time at coord even if it is not on a point of the grid
double getDistance(Coordinate coord1,Coordinate coord2);
// Returns the distance between two points
void calculateShiftForEachShot(VelocityModel* velModel, Configuration* config);
// The eikonal is far more precise when the source is set on a grid point. This function calculate the optimum xminGrid,yminGrid
// zminGrid in that purpose
std::vector<double> downSampleProfile(std::vector<double> profile, std::vector<double> zp, std::vector<double> zFiltp);
/* Downsample (or upsample) a curve defined as :
  |       v[0]       |       v[1]       |        ...        |      v[nz-1]      |
z[0]               z[1]               z[2]                z[nz-1]             z[nz]
           |                  |        ...        |                   |
           v                  v                   v                   v
         zp[0]              zp[1]              zp[nz-2]              zp[nz-1]

To obtain :
  |         vp[0]         |         vp[0]        |  ... |     vp[nzFilt-1]      |
zFilt[0]                zFilt[1]               zFilt[2] zFilt[nzFilt-1]      zFilt[nzFilt]
              |                       |     ...      |              |
              v                       v              v              v
           zFiltp[0]                 zFiltp[1]  zFiltp[nzFilt-2] zFiltp[nzFilt-1]

If nzFilt < nz => downsample 
If nzFilt > nz => upsample
 
vp in intrapolated by a linear approximation */
void InverseWaveletTransform(std::vector<double>* filteredLog, const std::vector<double>* coeffsToKeep, Configuration* config, int wavelet, bool sWaves);
// Calculate the P or S waves velocity profile corresponding to the wavelet coefficients "coeffsToKeep" (P waves : sWaves=false)
void InverseLayerTransform(std::vector<double>* filteredProfile, const std::vector<double>* params, Configuration* config, bool sWaves);
// Calculate the P or S waves velocity profile corresponding to the layers params (P waves : sWaves=false)
void meshing(std::vector<double>* profile, VelocityModel* velModel, bool swaves);
// Extend the velocity profile on the whole mesh (coarsen the sampling)
double energy(State* state,Chain* chain, Configuration* config);
// Give the energy of a state in a given chain (we need to give the chain to get the temperature)
void makeVel(VelocityModel* velModel, State* state, const Configuration* config);
// Build the velocity model corresponding to the parameters of the state
void eikonal3d(tab3d<double>* tt3d,const VelocityModel* velModel,const Configuration* config, int shotNumber, bool sWaves);
// Calculate the P or S waves (sWaves==false or sWaves==true) travel times corresponding to the shot number i the velocity model and
// the geometric configuration stored in config. Call the fortran subroutine : FTeik3d.f90 which needs Include_FTeik3d_2.0.f
std::vector<int> makeIndex(Chain* chain, Configuration* config);
// Make the index of the parameters to modify
// We could have done far better to deal with the static parameters but I was tired...
void priorIteration(Chain* chain,Configuration* config);
// Perform an iteration in the prior
void iterationMHindependent(Chain* chain,Configuration* config);
// Perform an independent MH iteration on the chain given (Often the highest temperature chain)
void iterationMH(Chain* chain,Configuration* config);
// Perform a classical MH iteration on the chain given.
void importanceSamplingSwap(Run* run, int i, Configuration* config);
// IR swap iteration, i is the index of the chain considered
int pickastate2(Run* run, int i, Configuration* config);
// Pick a state of the i th chain of the run according to an IS (importance sampling) draw
void updateSCI(Run* run, int n);
// Add a new line to SCI (importance weights (cumulative sum) and normalization coefficients)
void updateMinMaxProfiles(Chain* chain, Configuration* config);
// Update in and max velocities investigated
// void updateAverageProfiles(Chain* chain,Configuration* config, int i); old version
// Update the average, variance and quantiles profiles
void updateAverageProfiles(Run* run,Configuration* config, int i);
// Update the average, variance and quantiles profiles
void isThisConfigOk(VelocityModel* velModel, Configuration* config);
// Set refVelModel->OK to 1 if we are precise enough with the eikonal and this config
void finalizeRun(Run* run);
// Free memory allocated for the chains and finalize MPI

#endif /* FUNCTIONS_H_ */
