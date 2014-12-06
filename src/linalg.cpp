
/* File: linalg.cpp */

#include "linalg.h"

int fact(int n){

  if (n < 0)
    return 0;
  if (n == 0)
    return 1;
  else{
    return n * fact(n-1); 
  }
}


double rms(double* seq, int n){

  double avg = 0.0;

  for (int i = 0; i < n; i++){
    avg += seq[i];
  }

  avg = avg/n;
  
  double rms = 0.0;

  for (int i = 0; i < n; i++){
    rms += (seq[i] - avg) * (seq[i] - avg);
  }
  
  rms = sqrt(rms);

  seq[0] = 5.0;

  return rms;
}
