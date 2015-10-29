#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <errno.h>  // To detect mathematical errors
#include "wavelet2s.h"
#include <fftw3.h>
#include "functions.h"
#include "structures.h"
#include "defines.h"
#include "rngs.h"
#include "rvgs.h"
#include "generalFunctions.h"
#include "filesAndControl.h"

double getTime(tab3d<double>* tt3d,Coordinate coord,VelocityModel* velModel)
// Give the time at coord even if it is not on a point of the grid
{
  double Z=(coord.z-velModel->zmin)/velModel->dz +1.0; // Z=zin/dzin+1.
  double X=(coord.x-velModel->xmin)/velModel->dx +1.0; // X=xin/dxin+1.
  double Y=(coord.y-velModel->ymin)/velModel->dy +1.0; // Y=yin/dyin+1.

  int iz=(int)floor(Z), ix= (int)floor(X), iy= (int)floor(Y);
  int kz=iz+1, kx=ix+1, ky=iy+1;
  double dz1=Z-(double)iz,dx1=X-(double)ix,dy1=Y-(double)iy;
  double dz2=1.0-dz1,dx2=1.0-dx1,dy2=1.0-dy1;
  double t = dz2 * dx2 * dy2 * tt3d->get(iz-1,ix-1,iy-1)
          +  dz1 * dx2 * dy2 * tt3d->get(kz-1,ix-1,iy-1)
          +  dz2 * dx1 * dy2 * tt3d->get(iz-1,kx-1,iy-1)
          +  dz2 * dx2 * dy1 * tt3d->get(iz-1,ix-1,ky-1)
          +  dz2 * dx1 * dy1 * tt3d->get(iz-1,kx-1,ky-1)
          +  dz1 * dx2 * dy1 * tt3d->get(kz-1,ix-1,ky-1)
          +  dz1 * dx1 * dy2 * tt3d->get(kz-1,kx-1,iy-1)
          +  dz1 * dx1 * dy1 * tt3d->get(kz-1,kx-1,ky-1);
  return t;
}

std::vector<double> downSampleProfile(std::vector<double> v, std::vector<double> zp, std::vector<double> zFiltp)
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
 
vp in intrapolated by a linear approximation
!! warning !! the down sampled profile vFilt don't give the velocity at zFilt but between these points (at zFiltp). If zFilt has a size nFilt it 
will have a size nzFilt-1 (Likewise if z has a size nz, the given profile don't give the velocity at z but between these points : it has a 
size nz-1) */
{
  std::vector<double> vFilt;
  std::vector<int> idxInfSup(2,-1);
  int idxInf=-1,idxSup=-1;
  for (int i=0;i<(int)zFiltp.size();i++) { 
    idxInfSup=findInterval(zp,zFiltp[i]); // Returns the index idxInfSup={idxInf,idxSup} of vector verifying : 
                                          //   vector[idxInf] < value < vector[idxSup]
    idxInf=idxInfSup[0];
    idxSup=idxInfSup[1];
    vFilt.push_back((v[idxSup]-v[idxInf])/(zp[idxSup]-zp[idxInf])*(zFiltp[i]-zp[idxInf])+v[idxInf]);
  }  
  return vFilt;
}

void InverseWaveletTransform(std::vector<double>* filteredProfile, const std::vector<double>* params, Configuration* config, bool sWaves)
// Calculate the P or S waves velocity profile corresponding to the wavelet coefficients params (P waves : sWaves=false)
// Apparently config->flagP == config->flagS and config->lengthP == config->lengthS : hence we would not need the flag sWaves
// but we let it here just in case of.
{
  // ******** INVERSE DISCRETE WAVELET TRANSFORM IMPLEMENTATION*********"
  int numberOfCoeffsKept=0;
  if (sWaves == false) { // If we perform the wavelet transform to P waves velocity profile
    if(config->swaves) // If we are working with S waves as well
      numberOfCoeffsKept=(int)((double)config->data.indexParameters.size()/2.0); // If we are working with swaves, the first parameters are for P waves coefficients
    else // If we are just working with P waves
      numberOfCoeffsKept=(int)((double)config->data.indexParameters.size());
    for(int j=0;j<numberOfCoeffsKept;j++) // Loop on the coefficients
      config->coeffsP[(int)config->data.indexParameters[j]]=(*params)[j]; // Copy them into config->coeffsP
    idwt(config->coeffsP, config->flagP, config->wavelet, *filteredProfile, config->lengthP);  // Performs IDWT 
  }
  else { // If we perform the wavelet transform to S waves velocity profile
    if(!config->swaves) { // If we are not working with s waves
      std::cout << "Request for a wavelet transform of a S waves velocity profile while we should not. Terminating..." << std::endl;
      exit(0);
    }
    numberOfCoeffsKept=(int)((double)config->data.indexParameters.size()/2.0);
    for(int j=0;j<numberOfCoeffsKept;j++) // Loop on the coefficients
      config->coeffsS[(int)config->data.indexParameters[j]]=(*params)[j+numberOfCoeffsKept]; // Copy them into config->coeffsS
    idwt(config->coeffsS, config->flagS, config->wavelet, *filteredProfile, config->lengthS);  // Performs IDWT      
  }
}

void InverseLayerTransform(std::vector<double>* filteredProfile, const std::vector<double>* params, Configuration* config, bool sWaves)
// Calculate the P or S waves velocity profile corresponding to the layers params (P waves : sWaves=false)
{
  int numberOfCoeffsKept=0;
  if (sWaves == false) { // If we perform the layer transform to P waves velocity profile
    if(config->swaves) // If we are working with S waves as well
      numberOfCoeffsKept=(int)((double)config->data.minParameters.size()/2.0); // If we are working with swaves, the first parameters are for P waves coefficients
    else // If we are just working with P waves
      numberOfCoeffsKept=(int)((double)config->data.minParameters.size());
    for(int j=0;j<numberOfCoeffsKept;j++) // Loop on the coefficients
      config->coeffsP[j]=(*params)[j]; // Copy them into config->coeffsP
    // Store the filtered profile:
    (*filteredProfile).resize((int)config->data.zp.size());
    for(unsigned int i=0;i<config->data.zp.size();i++) {
      if ((int)i < config->indices[0])
        (*filteredProfile)[i] = config->coeffsP[0];
      else if ((int)i >= config->indices.back())
        (*filteredProfile)[i] = config->coeffsP.back();
      else {
        for (int j=0;j<config->npu-2;j++) {
          if ((int)i >= config->indices[j] and (int)i < config->indices[j+1])
            (*filteredProfile)[i] = config->coeffsP[j+1];
        }
      }
    }
  }
  else { // If we perform the layer transform to S waves velocity profile
    if(!config->swaves) { // If we are not working with s waves
      std::cout << "Request for a layer transform of a S waves velocity profile while we should not. Terminating..." << std::endl;
      exit(0);
    }
    numberOfCoeffsKept=(int)((double)config->data.minParameters.size()/2.0);
    for(int j=0;j<numberOfCoeffsKept;j++) // Loop on the coefficients
      config->coeffsS[j]=(*params)[j+numberOfCoeffsKept]; // Copy them into config->coeffsS
    // Store the filtered profile:
    (*filteredProfile).resize((int)config->data.zp.size());
    for(unsigned int i=0;i<config->data.zp.size();i++) {
      if ((int)i < config->indices[0])
        (*filteredProfile)[i] = config->coeffsS[0];
      else if ((int)i >= config->indices.back())
        (*filteredProfile)[i] = config->coeffsS.back();
      else {
        for (int j=0;j<config->npu-2;j++) {
          if ((int)i >= config->indices[j] and (int)i < config->indices[j+1])
            (*filteredProfile)[i] = config->coeffsS[j+1];
        }
      }
    }
  }
}

void meshing(std::vector<double>* profile, VelocityModel* velModel, bool swaves) 
// Extend the velocity profile on the whole mesh vel
{
  if((int)profile->size() != velModel->nz-1) {
    std::cout << "The profile can't be extended to the whole mesh... not the same number of z points! (" << (int)profile->size();
    std::cout << " != " << velModel->nz-1 << ")" << std::endl;
    exit(0);
  }
  for(int iz=0; iz < velModel->nz-1; iz++) {
    for(int ix=0; ix < velModel->nx-1; ix++) {
      for(int iy=0; iy < velModel->ny-1; iy++) {
        if(swaves)
          velModel->velS->set(iz,ix,iy,(*profile)[iz]);
        else
          velModel->velP->set(iz,ix,iy,(*profile)[iz]);
      }
    }
  }
}

void makeVel(VelocityModel* velModel, State* state, Configuration* config)
// Build the velocity model corresponding to the parameters of the state : this will have to be change each time we change the problem.
{
  std::vector<double> filteredProfileP;
  if (config->waveletParameterization)
    InverseWaveletTransform(&filteredProfileP,&(state->params), config, false); // Calculate the P waves velocity profile corresponding to the wavelet coefficients (nz-1 points)
  else
    InverseLayerTransform(&filteredProfileP,&(state->params), config, false); // Calculate the P waves velocity profile corresponding to the layers coefficients (nz-1 points)
  std::vector<double> downSampledPvel = downSampleProfile(filteredProfileP,config->data.zp,config->data.zFiltp); // downSampledPvel Contains the down sampled filtered velocity model !! It will have a size nzFilt-1
  for(int i=0; i < (int)downSampledPvel.size(); i++) {
    if(downSampledPvel[i]<=0.0) // A velocity can't be < 0 ! If it is the case...
      downSampledPvel[i] = 0.01; // We put the velocity to 0.01m/s on this point -> the travel times will be huge and the model will be rejected
  }
  meshing(&downSampledPvel,velModel,false); // Extend this profile on the whole mesh
  if(config->swaves) {
    std::vector<double> filteredProfileS;
    if (config->waveletParameterization)
      InverseWaveletTransform(&filteredProfileS,&(state->params), config, true); // Calculate the S waves velocity profile corresponding to the wavelet coefficients (nz-1 points)
    else
      InverseLayerTransform(&filteredProfileS,&(state->params), config, true); // Calculate the S waves velocity profile corresponding to the layers coefficients (nz-1 points)
    std::vector<double> downSampledSvel = downSampleProfile(filteredProfileS,config->data.zp,config->data.zFiltp); // downSampledSvel Contains the down sampled filtered velocity model !! It will have a size nzFilt-1
    for(int i=0; i < (int)downSampledSvel.size(); i++) {
      if(downSampledSvel[i]<=0.0) // A velocity can't be < 0 ! If it is the case...
        downSampledSvel[i] = 0.01; // We put the velocity to 0.01m/s on this point -> the travel times will be huge and the model will be rejected
    }
    meshing(&downSampledSvel,velModel,true); // Extend this profile on the whole mesh (coarsen the sampling)
  }
}

void eikonal3d(tab3d<double>* tt3d, const VelocityModel* velModel,const Configuration* config, int shotNumber, bool sWaves)
// Calculate the P or S waves (sWaves==false or sWaves==true) travel times corresponding to the shot number i the velocity model and
// the geometric configuration stored in config.
{
  int nx = tt3d->get_nx(), ny = tt3d->get_ny(), nz = tt3d->get_nz();
  float dzin = velModel->dz, dxin = velModel->dx, dyin = velModel->dy;
  int nsweep = config->nSweeps;
  float epsin = config->epsin ;
  float zsin=(float)config->data.coordShots[shotNumber].z-(float)velModel->zmin;
  float xsin=(float)config->data.coordShots[shotNumber].x-(float)velModel->xmin;
  float ysin=(float)config->data.coordShots[shotNumber].y-(float)velModel->ymin;
  //std::cout << "nx :" << nx << " ny : " << ny << " nz : " << nz << std::endl;
  //std::cout << "dxin :" << dxin << " dyin : " << dyin << " dzin : " << dzin << std::endl;
  //std::cout << "xsin :" << xsin << " ysin : " << ysin << " zsin : " << zsin << std::endl;
  //std::cout << "nsweep :" << nsweep << " epsin : " << epsin << std::endl;
  if(sWaves==true) 
    fteik_(velModel->velS->get_values(),tt3d->get_values(),&nz,&nx,&ny,&zsin,&xsin,&ysin,&dzin,&dxin,&dyin,&nsweep,&epsin); // Call to Fortran function
  else
    fteik_(velModel->velP->get_values(),tt3d->get_values(),&nz,&nx,&ny,&zsin,&xsin,&ysin,&dzin,&dxin,&dyin,&nsweep,&epsin); // Call to Fortran function
}

double energy(State* state,Chain* chain,Configuration* config)
// Give the energy of a state in a given chain (we need to give the chain to get the temperature)
// T0 is calculated just using P waves TODO
{
  double E=0.0;
  ArrivalTimes arrivalTimes;
  int numberOfShots = (int)config->data.coordShots.size();
  int numberOfStations = (int)config->data.coordStations.size();
  VelocityModel velModel;
  velModel.nx=config->data.nx; velModel.ny=config->data.ny; velModel.nz=config->data.nzFilt;
  velModel.dx=config->data.dx; velModel.dy=config->data.dy; velModel.dz=config->data.dzFilt;
  velModel.xmin = config->data.xminGrid; velModel.ymin = config->data.yminGrid; velModel.zmin = config->data.zminGrid;
  velModel.velP= new tab3d<double>(velModel.nz-1,velModel.nx-1,velModel.ny-1,-1);  // Create a (nz-1,nx-1,ny-1) mesh and initialize every cell at -1
  if (config->swaves)
    velModel.velS= new tab3d<double>(velModel.nz-1,velModel.nx-1,velModel.ny-1,-1);  // Create a (nz-1,nx-1,ny-1) mesh and initialize every cell at -1
  makeVel(&velModel,state,config); // Build the velocity model corresponding to the parameters of the state for P waves
  // We save the profile into chain->profilesP and chain->profilesS
  for(int iz=0;iz<config->data.nzFilt-1;iz++) {
    chain->profilesP[iz].push_back(velModel.velP->get(iz,1,1)); // Let us save the vertical line : (:,1,1) the velocity is 1D anyway
    if(config->swaves)
      chain->profilesS[iz].push_back(velModel.velS->get(iz,1,1)); // Let us save the vertical line : (:,1,1) the velocity is 1D anyway
  }
  if (chain->maxP.size() > 0)
    updateMinMaxProfiles(chain,config); // Update min and max velocities investigated

  if(config->test==1) // We don't enter the Eikonal, we generate random energies uniformly distributed between config->minEtest and config->maxEtest
    E=Uniform(config->minEtest,config->maxEtest);
  else {
    tab3d<double> tt3dP(config->data.nzFilt,config->data.nx,config->data.ny,-1.0);
    tab3d<double> tt3dS(config->data.nzFilt,config->data.nx,config->data.ny,-1.0); // TODO : this is declared even if swaves == false -> change that
    double sumP=0.0, sumS=0.0, totalSumP=0.0, totalSumS=0.0;
    std::vector<int> indexP, indexS;
    double t0=0.0;   // Origin time of the shot (if config->recalculateT0 == 1 we will recalculate it)
    int numberOfEikonalToCompute = 0;

    if(config->swaves) 
      numberOfEikonalToCompute = 2*numberOfShots;
    else
      numberOfEikonalToCompute = numberOfShots;
      
    #ifdef PAR
      MPI_Barrier(MPI_COMM_WORLD);
    #endif      
    
    for(int i=config->mpiConfig.rank;i<numberOfEikonalToCompute;i+=config->mpiConfig.nb_process) { // Loop on the shots
      if (i<numberOfShots) { // P waves
        if(config->verbose2)
          std::cout << "     Eikonal P for shot number " << i+1 << " on " << numberOfShots << " ..."<< std::endl;
        indexP.push_back(i);
 
        eikonal3d(&tt3dP,&velModel,config,i,false); 
 
        //Calculate the P waves travel times everywhere on the mesh (put them on tt3dP) for the shot number i
        for(int j=0;j<numberOfStations;j++) 
          arrivalTimes.timesP.push_back(getTime(&tt3dP,config->data.coordStations[j],&velModel));

        if (config->recalculateT0) {
          t0=0.0;
          int nUsedP=0;  // Number of travel times used for P waves
          for(int k=0;k<(int)arrivalTimes.timesP.size();k++) { // Calculation of t0
            int idxP=indexP[(int)floor(k/numberOfStations)]*numberOfStations+k-numberOfStations*(int)floor(k/numberOfStations); // index of the arrival time
            if (config->data.times.timesP[idxP] > 0.0) {
              t0+=config->data.times.timesP[idxP]-arrivalTimes.timesP[k];
              nUsedP++;
            }
          }
          t0=t0/nUsedP;
   //       #ifdef PAR
            // Send t0 to the process that will calculate the eikonal S for the same shot :
   //         MPI_Send(&t0, 1, MPI_DOUBLE, config->mpiConfig.rank+numberOfShots, config->mpiConfig.rank, MPI_COMM_WORLD);
   //       #endif
   //          std::cout << "     t0 for shot " << i+1 << " (P waves) : "<< t0 << std::endl;
        }
      } 
      else { // S waves
        if(config->verbose2)
          std::cout << "     Eikonal S for shot number " << i+1-numberOfShots << " on " << numberOfShots << " ..."<< std::endl;
        indexS.push_back(i-numberOfShots);
        eikonal3d(&tt3dS,&velModel,config,i-numberOfShots,true); 
        // Calculate the S waves travel times everywhere on the mesh (put them on tt3dS) for the shot number i
        for(int j=0;j<numberOfStations;j++) 
          arrivalTimes.timesS.push_back(getTime(&tt3dS,config->data.coordStations[j],&velModel));
       
        if (config->recalculateT0) {
          t0=0.0;
          int nUsedS=0;  // Number of travel times used for P waves
          for(int k=0;k<(int)arrivalTimes.timesS.size();k++) { // Calculation of t0
            int idxS=indexS[(int)floor(k/numberOfStations)]*numberOfStations+k-numberOfStations*(int)floor(k/numberOfStations); // index of the arrival time
            if (config->data.times.timesS[idxS] > 0.0) {
              t0+=config->data.times.timesS[idxS]-arrivalTimes.timesS[k];
              nUsedS++;
            }
          }
          t0=t0/nUsedS;
   //       std::cout << "     t0 for shot " << i+1-numberOfShots << " (S waves) : " << t0 << std::endl;
        }
    //    #ifdef PAR
    //      if (config->recalculateT0) {
    //        MPI_Status status;
            // Receive t0 from the process that has calculated the eikonal P for the same shot :
    //        MPI_Recv(&t0, 1, MPI_DOUBLE,  config->mpiConfig.rank-numberOfShots, config->mpiConfig.rank-numberOfShots, MPI_COMM_WORLD, &status);
    //      }
    //    #endif
      }  
    }

    for(int k=0;k<(int)arrivalTimes.timesP.size();k++) { // Loop on the P wave travel times calculated
      int idxP=indexP[(int)floor(k/numberOfStations)]*numberOfStations+k-numberOfStations*(int)floor(k/numberOfStations); // index of the arrival time corresponding to the one calculated. 
      double diffP=0.0;
      if (config->data.times.timesP[idxP] > 0.0)
        diffP= config->data.times.timesP[idxP]-(t0+arrivalTimes.timesP[k]);
   //   std::cout << " config->data.times.timesP[" << idxP << "]  : " << config->data.times.timesP[idxP] << "  arrivalTimes.timesP[" << k << "] : " << arrivalTimes.timesP[k] << "  Diff : " << diffP << std::endl; 
      sumP+=pow(diffP/config->data.sigmaP,2.0)/2.0; // sum of the squares of the differences between the P wave first arrival times and the data
    }
    if(config->swaves) { 
      for(int k=0;k<(int)arrivalTimes.timesS.size();k++) { // Loop on the S wave travel times calculated
        int idxS=indexS[(int)floor(k/numberOfStations)]*numberOfStations+k-numberOfStations*(int)floor(k/numberOfStations); // index of the arrival time corresponding to the one calculated
        double diffS=0.0;
        if (config->data.times.timesS[idxS] > 0.0)
          diffS= config->data.times.timesS[idxS]-(t0+arrivalTimes.timesS[k]);
    //    std::cout << " config->data.times.timesS[" << idxS << "]  : " << config->data.times.timesS[idxS] << "  arrivalTimes.timesS[" << k << "] : " << arrivalTimes.timesS[k] << "  Diff : " << diffS << std::endl;   
        sumS+=pow(diffS/config->data.sigmaS,2.0)/2.0; // sum of the squares of the differences between the P wave first arrival times
      }
    }
    #ifdef PAR
      MPI_Barrier(MPI_COMM_WORLD);
      if(config->swaves) { 
        MPI_Allreduce(&sumP,&totalSumP,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD) ; 
        MPI_Allreduce(&sumS,&totalSumS,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD) ; 
      }
      else 
        MPI_Allreduce(&sumP,&totalSumP,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD) ; 
    #else 
      totalSumP=sumP;
      totalSumS=sumS;
    #endif
    E=(totalSumP+totalSumS+config->data.Ep)/(chain->T);
    if (config->verbose1 && config->mpiConfig.rank == 0) { 
      std::cout << "      Energy : " << E << "  (sumP : " << totalSumP << "  sumS : " << totalSumS << "  Ep : " << config->data.Ep << ")  E*T : " << E*chain->T << ")" << std::endl;
    }
  }
  delete velModel.velP;
  if (config->swaves)
    delete velModel.velS;
  return E;
  getchar();
}

std::vector<int> makeIndex(Chain* chain, Configuration* config)
// Make the index of the parameters to modify
// We could have done far better to deal with the static parameters but I was tired...
{
  int nParams=(int)chain->states.back().params.size(); // = NPU (or NPU*2)
  int varyingParams=nParams; // Number of parameters that will actually vary
  std::vector<int> index;
  for(int i=0;i<nParams;i++)
    index.push_back(i);
  for (int i=0;i<(int)config->data.staticParameters.size();i++) { // Remove the static parameters from the indexes that could be modified
    index.erase(index.begin()+config->data.staticParameters[i]-i);   // Erase the static elements
    varyingParams--;
  }
  for(int j=varyingParams;j>chain->nc;j--) { // At the end of the loop index contains a random series of nc indexes to modify
    double r=((j-0.5)-(0+0.5))*rand()/(double)RAND_MAX+1;  // To draw randomly the parameter to modify : c
    int c = (int)(r < 0 ? r - 0.5 : r + 0.5);
      index.erase(index.begin()+c-1);   // Erase the c-1 element
  }
  if(config->verbose2 && config->mpiConfig.rank == 0)
  {
    std::cout << std::endl << "Parameters that will be modified : " << std::endl;
    for(int k=0;k<(int)index.size();k++)
      std::cout << index[k] << " ";
    std::cout << std::endl;
  }
  return index;
}

void priorIteration(Chain* chain,Configuration* config)
// Perform an iteration in the prior
{
  if (config->verbose1 && config->mpiConfig.rank == 0) {
    std::cout << "***** Iteration on chain " << chain->i+1 << " on " << config->nbt << " (Temperature : "<< config->T[chain->i] << ") ***** " << std::endl;
    std::cout<< "--> Prior iteration... " << std::endl;
  }
  State priorState=chain->states.back();
  // Copy the actual state into the trial one. chain->states.back() is the actual state
  std::vector<int> index = makeIndex(chain,config); // Draw randomly which parameters will be modified
  for(int k=0;k<(int)index.size();k++) {  // Loop on the number of parameters that will be modified
    priorState.params[index[k]]=Uniform(config->data.minParameters[index[k]],config->data.maxParameters[index[k]]);
    if (config->verbose2 && config->mpiConfig.rank == 0) {
      std::cout<< "priorState.params["<< index[k] <<"] : " << priorState.params[index[k]] << "   ";
      std::cout<< "(Uniform between : "<< config->data.minParameters[index[k]] <<" and " << config->data.maxParameters[index[k]] << ") " << std::endl;
    }
  }
  priorState.E=energy(&priorState,chain,config); // Compute the forward problem for the trial model
  chain->states.push_back(priorState);
}

void iterationMHindependent(Chain* chain,Configuration* config)
// Perform an independent MH iteration on the chain given (Often the highest temperature chain)
// We could have done far better to deal with the static parameters but I was tired... TODO : add a function
{
  if (config->verbose1 && config->mpiConfig.rank == 0) {
    std::cout << "***** Iteration on chain " << chain->i+1 << " on " << config->nbt << " (Temperature : "<< config->T[chain->i] << ") ***** " << std::endl;
    std::cout<< "--> Independent Metropolis Hasting iteration proposed... " << std::endl;
  } 
  State trialState=chain->states.back();
  // Copy the actual state into the trial one. chain->states.back() is the actual state
  std::vector <int> index;
  index=makeIndex(chain,config); // Draw randomly which parameters will be modified
  for(int k=0;k<(int)index.size();k++) { // Loop on the number of parameters that we will modify
    trialState.params[index[k]]=Uniform(config->data.minParameters[index[k]],config->data.maxParameters[index[k]]);
    if (config->verbose2 && config->mpiConfig.rank == 0 && 0) {
      std::cout<< "Trial state param["<< index[k] <<"] : " << trialState.params[index[k]] << "   ";
      std::cout<< "(Uniform between : "<< config->data.minParameters[index[k]] <<" and " << config->data.maxParameters[index[k]] << ") " << std::endl;
    }  
  }
  trialState.E=energy(&trialState,chain,config); // Compute the forward problem for the trial model (add one line to the profiles)
  if (config->verbose2 && config->mpiConfig.rank == 0) {
    std::cout<< "trialState energy :" << trialState.E << std::endl;
    std::cout<< "Last state energy :" << chain->states.back().E << std::endl;
    std::cout<< "alpha :" << exp(chain->states.back().E-trialState.E) << std::endl;
  }
  double p=Random();  // Acceptance variable
  long double alpha=exp(chain->states.back().E-trialState.E);
  // Acceptance probability alpha=exp(E[n-1][i]-Ep)
  if (p<alpha)  // We accept the transition
  {
    if (config->verbose1 && config->mpiConfig.rank == 0)
      std::cout<< "Transition accepted" << std::endl<< std::endl;
    chain->states.push_back(trialState);
    chain->at++; // Number of accepted transitions
  }
  else // We keep the same model
  {
    if (config->verbose1 && config->mpiConfig.rank == 0)
      std::cout<< "Transition rejected" << std::endl<< std::endl; 
    chain->states.push_back(chain->states.back());
    chain->rt++; // Number of rejected transitions
    for(int iz=0;iz<config->data.nzFilt-1;iz++) {
      chain->profilesP[iz].pop_back(); //.erase(chain->profilesP[iz].end() - 1); // Erase last element
      chain->profilesP[iz].push_back(chain->profilesP[iz].back());
      if(config->swaves) {
        chain->profilesS[iz].pop_back(); //.erase(chain->profilesS[iz].end() - 1); // Erase last element
        chain->profilesS[iz].push_back(chain->profilesS[iz].back());
      }
    }
  }
}

void iterationMH(Chain* chain,Configuration* config)
// Perform a classical MH iteration on the chain given.
// We could have done far better to deal with the static parameters but I was tired... TODO : add a function
{
  if (config->verbose1 && config->mpiConfig.rank == 0) {
    std::cout << "***** Iteration on chain " << chain->i+1 << " on " << config->nbt << " (Temperature : "<< config->T[chain->i] << ") ***** " << std::endl;
    std::cout<< "--> Metropolis Hasting iteration proposed... " << std::endl;
  }
  State trialState=chain->states.back();
  // Copy the actual state into the trial one. chain->states.back() is the actual state
  std::vector <int> index;
  index=makeIndex(chain,config); // Make the index of the parameters to modify
  bool outOfDomain = false;
  double sigma = 0.0;
  double oldValue = 0.0;
  for(int k=0;k<(int)index.size();k++) { // Loop on the number of parameters that we will modify
    ///////////IMPORTANT LINES////////////
    oldValue = chain->states.back().params[index[k]];
    // double sigma = chain->deltaParameters[index[k]]*sqrt(chain->T/config->tmax);
    // double sigma = chain->deltaParameters[index[k]]*pow((double)(chain->i+1)/(double)config->nbt,config->n);
    // double sigma = chain->deltaParameters[index[k]]*sqrt(chain->T/config->tmax)*floor((double)config->nbt/(double)(chain->i+1));
    sigma = chain->deltaParameters[index[k]]; 
    trialState.params[index[k]]=oldValue+Normal(0.0,sigma);
    /////////////////////////////////////
    if (config->verbose2 && config->mpiConfig.rank == 0) {
      std::cout<< "  Trial state param["<< index[k] <<"] : " << trialState.params[index[k]] << "   ";
      std::cout<< "  Last state param["<< index[k] <<"] : " << chain->states.back().params[index[k]] << "   ";
      std::cout<< "  (Gaussian centred in : " << oldValue << " and of delta : " << sigma << " )" << std::endl;
      std::cout<< "  Need to be between : "<< config->data.minParameters[index[k]] <<" and " << config->data.maxParameters[index[k]] << std::endl;
    } 
    if (trialState.params[index[k]] > config->data.maxParameters[index[k]] || trialState.params[index[k]] < config->data.minParameters[index[k]])
      outOfDomain=true;
  }
  if (outOfDomain==false) { // We can iterate
    trialState.E=energy(&trialState,chain,config); // Compute the forward problem for the trial model (add one line to the profiles)
    if (config->verbose2 && config->mpiConfig.rank == 0) {
      std::cout<< "trialState energy :" << trialState.E << std::endl;
      std::cout<< "Last state energy :" << chain->states.back().E << std::endl;
      std::cout<< "alpha :" << exp(chain->states.back().E-trialState.E) << std::endl;
    }
    double p=Random();  //Acceptance variable
    long double alpha=exp(chain->states.back().E-trialState.E);
    // Acceptance probability alpha=exp(E[n-1][i]-Ep)
    if (p<alpha) {  // We accept the transition
      if (config->verbose1 && config->mpiConfig.rank == 0)
         std::cout<< "Transition accepted" << std::endl<< std::endl;
      chain->states.push_back(trialState);
      chain->at++; // Number of accepted transitions
    }
    else { // We keep the same model
      if (config->verbose1 && config->mpiConfig.rank == 0)
        std::cout<< "Transition rejected" << std::endl<< std::endl;
      chain->states.push_back(chain->states.back());
      chain->rt++; // Number of rejected transitions
      for(int iz=0;iz<config->data.nzFilt-1;iz++) {
        chain->profilesP[iz].pop_back(); //.erase(chain->profilesP[iz].end() - 1); // Erase last element
        chain->profilesP[iz].push_back(chain->profilesP[iz].back());
        if(config->swaves) {
          chain->profilesS[iz].pop_back(); //.erase(chain->profilesS[iz].end() - 1); // Erase last element
          chain->profilesS[iz].push_back(chain->profilesS[iz].back());
        }
      }
    }
  }
  else { // Out of domain. We stay in the same state
    if (config->verbose1 && config->mpiConfig.rank == 0)
      std::cout<< "Out Of Domain" << std::endl<< std::endl;
    chain->states.push_back(chain->states.back());
    chain->od++; // Number of "out of domain"
    for(int iz=0;iz<config->data.nzFilt-1;iz++) {
      chain->profilesP[iz].push_back(chain->profilesP[iz].back());
      chain->profilesS[iz].push_back(chain->profilesS[iz].back());
    }
   /*
      chain->profilesP[iz].pop_back(); //.erase(chain->profilesP[iz].end() - 1); // Erase last element
      if(config->swaves) {
        chain->profilesS[iz].pop_back(); //.erase(chain->profilesS[iz].end() - 1); // Erase last element
      }
    }
    */
  }
}

void importanceSamplingSwap(Run* run,int i,Configuration* config)
// IR swap iteration, i is the index of the chain considered
{

  if (config->verbose1 && config->mpiConfig.rank == 0) {
    std::cout << "***** Iteration on chain " << i+1 << " on " << config->nbt << " (Temperature : "<< config->T[i] << ") ***** " << std::endl;
    std::cout<< "--> Importance-Sampling swap proposed... " << std::endl;
  //  std::cout<< "Energies of the states of the higher temperature chain : " << std::endl;
  //  for(int j=0;j<(int)run->chains[i+1]->states.size();j++) {
  //    std::cout << "Ei+1(" << j << ") : " << run->chains[i+1]->states[j].E << std::endl;
  //  }
  }
  // We choose a higher temperature chain's state according to an IS (importance sampling) draw
  run->chains[i]->ps++;  // Swapping suggested
  int sp=pickastate2(run,i,config); // Pick a state of the i th chain of the run according to an IS (importance sampling) draw
  int n=(int)run->chains[i]->states.size()-1; // Index of the last state of the chain i
  double dT= 1/run->chains[i]->T-1/run->chains[i+1]->T;
  double dE=run->chains[i]->states.back().E*(run->chains[i]->T)-run->chains[i+1]->states[sp].E*(run->chains[i+1]->T);  // ->states[n-1].E!
  long double alpha=exp(dT*dE);
  //alpha=exp((1/T[i+1]-1/T[i])*(E[sp][i+1]T[i+1]-E[n-1][i]T[i]));
  if (config->test)
    alpha=0.5;
  if (config->verbose2 && config->mpiConfig.rank == 0) {
    std::cout<< "There are " << run->chains[i+1]->states.size() << "  states on the higher temperature chain " << std::endl;
    std::cout<< "The state nb : " << sp << "  (energy : " << run->chains[i+1]->states[sp].E << " ) in chain "<< i+2 <<" has been picked" << std::endl;
    std::cout<< "Energy of the state picked on chain "<< i+2 << " * Ti+1/Ti : " << run->chains[i+1]->states[sp].E*run->chains[i+1]->T/run->chains[i]->T << std::endl;
    std::cout<< "Actual energy : " << run->chains[i]->states.back().E << std::endl;
    std::cout<< "alpha :" << alpha << std::endl;
  }
  double p=Random();
  if (p<alpha) { // We swap the states and we iterate with a classical MH
    if (config->verbose1 && config->mpiConfig.rank == 0)
      std::cout<< "Swapping accepted" << std::endl<< std::endl;
    run->chains[i]->as++;
    run->chains[i]->states.push_back(run->chains[i+1]->states[sp]);  
    // !! The two states are not at the same temperature so the energy need to be corrected !!
    run->chains[i]->states.back().E=run->chains[i]->states.back().E*run->chains[i+1]->T/run->chains[i]->T;
    for(int iz=0;iz<config->data.nzFilt-1;iz++) {
      run->chains[i]->profilesP[iz].push_back(run->chains[i+1]->profilesP[iz][sp]);
      if(config->swaves) {
        run->chains[i]->profilesS[iz].push_back(run->chains[i+1]->profilesS[iz][sp]);
      }
    }
// TODO update profiles
    /* TODO
    if (config->verbose1 && config->mpiConfig.rank == 0)
      std::cout<< "Iteration of the swapped state... " << std::endl;
    iterationMH(run->chains[i],config);  // We iterate the chain
    // We erase the temporary state
    run->chains[i]->states.erase(run->chains[i]->states.end()-1); // Delete the penultimate state
    */
    SwapFeatures swapFeatures;
    swapFeatures.idxChain=i; // Index of the chain which has swapped (0 : lowest temperature chain...)
    swapFeatures.idxState=n; // Index of the state of the chain who has been swapped
    swapFeatures.sp=sp;      // Index of the state of the higher temperature chain who has been chosen for the swap
    run->swapHist.push_back(swapFeatures); // Add the features of the swap to the swapping history off the run.
    writeSwap(run,config);
    // Writes a new line on the data file ll.XXX.dat. If it does not exist it is created. It will store the swapping history
  }
  else { // We keep the previous state
    if (config->verbose1 && config->mpiConfig.rank == 0)
      std::cout<< "Swapping rejected" << std::endl<< std::endl;
    run->chains[i]->states.push_back(run->chains[i]->states.back());
    run->chains[i]->rs++;
    for(int iz=0;iz<config->data.nzFilt-1;iz++) {
      run->chains[i]->profilesP[iz].push_back(run->chains[i]->profilesP[iz].back());
      if(config->swaves) {
        run->chains[i]->profilesS[iz].push_back(run->chains[i]->profilesS[iz].back());
      }
    }
  }
}

int pickastate2(Run* run, int i, Configuration* config)
// Pick a state of the i th chain of the run according to an IS (importance sampling) draw
{
  int r=0;
  int n=(int)run->chains[i]->states.size();
  long double p=(long double)Random()*(run->SCI[n-1][i]);

  // Draw a number between 0 and the importance weight of the i th chain's previous state  (p=Random()*SCI[n-1][l]);
//  if (config->verbose1 && config->mpiConfig.rank == 0) {
//    std::cout<< "run->SCI["<<n-1<<"][" << i<<"] : " << run->SCI[n-1][i] << std::endl<< std::endl;
//    std::cout<< "p : "<< p << std::endl<< std::endl;
//  }

  for (int j=0;j<n;j++) {
    if (p >run->SCI[j][i]) {
//      if (config->verbose1 && config->mpiConfig.rank == 0) 
//        std::cout<< "p > run->SCI["<< j << "][" << i<< "]    (" << run->SCI[j][i] << ")  ->  r++" << std::endl;
      r++;
    }
    else {
//      if (config->verbose1 && config->mpiConfig.rank == 0) 
//        std::cout<< "p < run->SCI["<< j << "][" << i<< "]    (" << run->SCI[j][i] << ")  " << std::endl;
    }
  }
//  if (config->verbose1 && config->mpiConfig.rank == 0) 
//   std::cout<< "r to be returned : "<< r << std::endl;
  return r;
}

void updateSCI(Run* run, int n)
// Add a new line to SCI (importance weights (cumulative sum) and normalization coefficients)
{
  std::vector<long double> line(run->chains.size()-1);
  for (unsigned int j=0;j<run->chains.size()-1;j++) 
    line[j]=run->SCI[n][j]+(long double)expl(run->chains[j+1]->states[n+1].E*(1-run->chains[j+1]->T/run->chains[j]->T));
  //  line[j]=run->SCI[n][j]+exp(run->chains[j+1]->states[n+1].E*(1.-run->chains[j+1]->T/run->chains[j]->T));
  //SCI[n][j]=SCI[n-1][j]+exp(E[n][j+1]T[j+1]*(1/T[j+1]-1/T[j]));
  run->SCI.push_back(line);
}

void updateMinMaxProfiles(Chain* chain, Configuration* config)
// Update in and max velocities investigated
{
  double lastPvel = 0.0, lastSvel=0.0;

  for (int iz=0;iz<config->data.nzFilt-1;iz++) { // Loop on the values
    lastPvel = chain->profilesP[iz].back();
 //   if(iz == 5 && chain->i == 2)
 //     std::cout << "last Pvel : " << lastPvel << std::endl;
    if (lastPvel > chain->maxP[iz]) // If the last element is the biggest one
      chain->maxP[iz] = lastPvel;
    if (lastPvel < chain->minP[iz]) // If the last element is the smallest one
      chain->minP[iz] = lastPvel;
    if(config->swaves) {
      lastSvel = chain->profilesS[iz].back();
      if (lastSvel > chain->maxS[iz]) // If the last element is the biggest one
        chain->maxS[iz] = lastSvel;
      if (lastSvel < chain->minS[iz]) // If the last element is the smallest one
        chain->minS[iz] = lastSvel;
      }
  }
//  if(chain->i == 2) {
//    std::cout << "maxP 5 for chain 2 : " << chain->maxP[5] << std::endl; // TODO remove
//    std::cout << "minP 5 for chain 2 : " << chain->minP[5] << std::endl; // TODO remove  
//  }
}

void updateAverageProfiles(Chain* chain,Configuration* config, int i)
// Update the average, variance and quantiles profiles
// TODO : put if config->swaves
{

/* The variance and the average are calculated with the algorithm :
average=0;
var=0;
for i=1:length(dataset)
   AverageOld=average;
   average=average+(dataset(i)-average)/i;
   var=var+(dataset(i)-averageOld)*(dataset(i)-average);
end
var=var/(length(dataset)-1);
*/
  int nPoints,nOut,idxInf,idxSup;
  std::vector<double> averageOldP, averageOldS;
  std::vector<double> sortedProfilesP, sortedProfilesS;
  for (int iz=0;iz<config->data.nzFilt-1;iz++) {
    averageOldP.push_back(chain->averageP[iz]);
    averageOldS.push_back(chain->averageS[iz]);
  }
// chain->profilesP[iz].back() -> Current P waves velocity at depth z
  for (int iz=0;iz<config->data.nzFilt-1;iz++) {
    chain->averageP[iz]=chain->averageP[iz]+(chain->profilesP[iz].back()-chain->averageP[iz])/((double)i+1.0);
    chain->averageS[iz]=chain->averageS[iz]+(chain->profilesS[iz].back()-chain->averageS[iz])/((double)i+1.0);
    chain->varP[iz]=chain->varP[iz]+(chain->profilesP[iz].back()-averageOldP[iz])*(chain->profilesP[iz].back()-chain->averageP[iz]);
    chain->varS[iz]=chain->varS[iz]+(chain->profilesS[iz].back()-averageOldS[iz])*(chain->profilesS[iz].back()-chain->averageS[iz]);
    nPoints=(int)chain->profilesP[iz].size(); // Number of profiles generated
    nOut=nPoints*(1-config->qp); // Number of points out of the quantile
    idxInf=floor(nOut/2); // index of the lowest velocity on the quantile
    idxSup=nPoints-ceil(nOut/2); // index of the highest velocity on the quantile
    sortedProfilesP=chain->profilesP[iz];
    sortedProfilesS=chain->profilesS[iz];
    std::sort (sortedProfilesP.begin(), sortedProfilesP.end());  // We sort the values (after having used the last value!)
    std::sort (sortedProfilesS.begin(), sortedProfilesS.end());  // We sort the values
    chain->qInfP[iz]=sortedProfilesP[idxInf];
    chain->qSupP[iz]=sortedProfilesP[idxSup];
    chain->qInfS[iz]=sortedProfilesS[idxInf];
    chain->qSupS[iz]=sortedProfilesS[idxSup];
  }

// chain->varP[iz]=chain->varP[iz]/((int)chain->ProfileP[iz].size()-1); -> Done when we write on files
// chain->varS[iz]=chain->varP[iz]/((int)chain->ProfileS[iz].size()-1); -> Done when we write on files

}

void finalizeRun(Run* run)
// Free memory allocated for the chains and finalize MPI
{
  for(unsigned int i=0;i<run->chains.size();i++)
    delete(run->chains[i]);          // Free the memory allocated for the chains
  #ifdef PAR
  MPI_Finalize();
  #endif
}

