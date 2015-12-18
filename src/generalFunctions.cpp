#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <sstream>
#include <functional>   // std::greater
#include <algorithm>    // std::sort, std::min_element, std::max_element
#include <stdlib.h>
#include <sys/stat.h>

#include "defines.h"
#include "structures.h"

#define SPACES " \t\r\n"

void convertTime(char * buffer, time_t rawtime)
// Convert a raw unix time into for example 03/Mar/2013 01:32:12. Store the string into a given char*
{
  struct tm * timeinfo;
  timeinfo = localtime ( &rawtime );
  strftime (buffer,80,"%d/%b/%Y %H:%M:%S.",timeinfo);
}

void printTime(Configuration* config)
{
  if (config->mpiConfig.rank == 0) {   // (In case of parallel implementation just one process has to create files)
    time_t t=time(NULL);
    if (t == -1)  // Sometime time() return -1
      std::cout << "TIME PROBLEM"<< std::endl;
    char buffer[80];
    convertTime(buffer, t);
    std::cout << "time : "<< buffer << std::endl;
  }
}

std::string formatBytes(unsigned long bytes)
// format bytes as a human readable string without branching
{
  char tmp[128] = "";
  const char *si_prefix[] = { "B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB" };
  const int base = 1024;
  int c = std::min((int)(log((double)bytes)/log((double)base)), (int)sizeof(si_prefix) - 1);
  sprintf(tmp, "%1.2f %s", bytes / pow((double)base, c), si_prefix[c]);
  return std::string(tmp);
}

std::vector<double> linspace(double first, double last, int len)
// Return an evenly spaced 1-d vector of doubles
{
  std::vector<double> result(len);
  double step = (last-first) / (len - 1);
  for (int i=0; i<len; i++) 
    result[i] = first + i*step;
  return result;
}

std::vector<int> findInterval(std::vector<double> vector, double value)
// Returns the index idxInfSup={idxInf,idxSup} of vector verifying : vector[idxInf] < value < vector[idxSup]
// If it is not possible returns -1,-1
{
  std::vector<int> idxInfSup(2,-1);
  std::vector<double> n;
  for(int i=0;i<(int)vector.size();i++)
    n.push_back(fabs(vector[i]-value));
  double minElem = *std::min_element(n.begin(), n.end()); // min distance to value
  int idx = std::find(n.begin(), n.end(), minElem) - n.begin(); // index of element that is closest to value
  double minVector = *std::min_element(vector.begin(), vector.end());
  double maxVector = *std::max_element(vector.begin(), vector.end());
  if ((value < minVector) || (value > maxVector)) {
    idxInfSup[0] = -1;
    idxInfSup[1] = -1;
  }
  else if ((vector[idx] >= value) && (idx != 0)) {
    idxInfSup[0] = idx-1;
    idxInfSup[1] = idx;
  }
  else {
    idxInfSup[0] = idx;
    idxInfSup[1] = idx+1;
  }
  return idxInfSup;
}

int sign(double x)
// If x > 0 sign(x) = 1, if x < 0 sign(x) = -1, if x = 0 sign(x) = 0
{
	return ((x < 0.0) ? -1 : ((x > 0.0) ? 1 : 0) );
}

void findthresh(std::vector<double> vector, int N, double& t)
// At the end t is the Nth biggest value in vector
{
  std::sort(vector.begin(), vector.end(), std::greater<double>());
  t = vector.at(N-1);
}

void keepNsignificantValues(std::vector<double>* vector, int N)
// Keep the N biggest absolute values in vector. Put the others to 0
{
  int notZeros = 0;
  double thresh = 0.0;
  std::vector<double> temp_dwtoutput;
  for (unsigned int i = 0; i < vector->size();i++) {
    double temp = fabs((*vector)[i]);
    temp_dwtoutput.push_back(temp);
  }
  findthresh(temp_dwtoutput,N, thresh); // thresh is the Nth biggest value of temp_dwtoutput
  for (unsigned int i = 0; i < vector->size();i++) {
    double temp = fabs((*vector)[i]);
    if (temp < thresh or notZeros >= N)
      (*vector)[i] = 0.0;
    else
      notZeros++;
  }
}

void keepNfirstValues(std::vector<double>* vector, int N)
// Keep the N first values in vector. Put the others to 0
// For example if vector = 10 21 12 12 7 8 and N = 2
// at the end vector = 10 21 12 0 0 0
{
  for (unsigned int i = 0; i < vector->size();i++) {
    if ((int)i > N-1)
      (*vector)[i] = 0.0;
  }
}

//*****Trim functions created by Fabien Ors*****
std::string trim_right (const std::string & s, const std::string & t = SPACES)
{
  std::string d (s);
  std::string::size_type i (d.find_last_not_of (t));
  if (i == std::string::npos)
    return "";
  else
    return d.erase (d.find_last_not_of (t) + 1) ;
}

std::string trim_left (const std::string & s, const std::string & t = SPACES)
{
  std::string d (s);
  return d.erase (0, s.find_first_not_of (t)) ;
}

std::string trim (const std::string & s, const std::string & t = SPACES)
{
  std::string d (s);
  return trim_left (trim_right (d, t), t) ;
}

int closest(std::vector<double>& vec, double value)
// Returns the indexes of the closest vector's value that's less than or equal to a given value
{
  std::vector<double>::iterator it = std::lower_bound(vec.begin(), vec.end(), value);
  if (it == vec.end()) { return -1; }
  return std::distance(vec.begin(), it);
}

const char* convertDoubleToChar(double value)
// From a double (ex:4561.54) return a char* (ex: "4561.54")
{
    std::ostringstream strs;
    strs.precision(PREC);
    strs << value;
    std::string strstr = strs.str();
    const char* valueChar = strstr.c_str();
    return valueChar;
}
