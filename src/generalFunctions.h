/*
 * generalFunctions.h
 *
 *  Created on: 13 may 2013
 *      Author: abottero
 */

#ifndef GENERALFUNCTIONS_H_
#define GENERALFUNCTIONS_H_
#define SPACES " \t\r\n"

#include "structures.h"

void convertTime(char * buffer, time_t rawtime);
// Convert a raw unix time into for example 12/Mar/2013 01:32:12. Store the string into a given buffer
void printTime(Configuration* config);
// Print the time. For example 12/Mar/2013 01:32:12
std::string formatBytes(unsigned long bytes);
// format bytes as a human readable string without branching
std::vector<double> linspace(double first, double last, int len);
// Return an evenly spaced 1-d vector of doubles
std::vector<int> findInterval(std::vector<double> vector, double value);
// Returns the index idxInfSup={idxInf,idxSup} of vector verifying : vector[idxInf] < value < vector[idxSup]
// If it is not possible returns -1,-1
int sign(double x);
// If x > 0 sign(x) = 1, if x < 0 sign(x) = -1, if x = 0 sign(x) = 0
void findthresh(std::vector<double> vector, int N, double& t);
// At the end t is the Nth biggest value in vector
void keepNsignificantValues(std::vector<double>* vector, int N);
// Keep the N biggest absolute values in vector. Put the others to 0
// !! Warning !! for example if vector = 10 21 12 12 7 8 and N = 2
// at the end vector = 0 21 12 12 0 0. We could fix that.

//*****Trim functions created by Fabien Ors*****
std::string trim_right (const std::string & s, const std::string & t = SPACES);
std::string trim_left (const std::string & s, const std::string & t = SPACES);
std::string trim (const std::string & s, const std::string & t = SPACES);
//**********************************************

#endif /* GENERALFUNCTIONS_H_ */
