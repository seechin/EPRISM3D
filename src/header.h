//#define HAVE_PTHREAD_H
//#define HAVE_ZLIB_H
#include "../config.h"

// Linux libraries
#ifdef HAVE_PTHREAD_H
  #define   _LOCALPARALLEL_                 // general local paralleling
  #define   _LOCALPARALLEL_PTHREAD_
#else
  #define   _LOCALPARALLEL_
#endif
#ifdef HAVE_ZLIB_H
  #define   _LIBZ_                          // use libz to compress TENSOR4D
#endif

// Third party libraries
// Gromacs for XTC
#ifdef GROMACS4
  #define   _GROMACS4_
#endif
#ifdef GROMACS16
  #define   _GROMACS2016_
#endif

// Features
//#define     _INTERACTIVE_                 // allow interactive mode, undef this will disable main-interactive.cpp
#define     _TTYPROMPTCOLOR_                // allow color text in terminal
//#define     _FUNCTION_EXPORT_
//#define     _FFTWMPPARALLEL_              // fftw multithreading
