#ifndef ANALYSIS_LIB_H
#define ANALYSIS_LIB_H

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>

extern char ALIB_line_identifier;


#define ALIB_INITIALIZATION_MODE 0
#define ALIB_SAMPLING_MODE 1
#define ALIB_PRINT_MODE 2
//#define ALIB_RETURN_DATA_MODE 3

#define ALIB_SUCCESS 0
#define ALIB_ERROR_NONE ALIB_SUCCESS
#define ALIB_ERROR  1
#define ALIB_ERROR_NO_DATA 2
#define ALIB_ERROR_NULL_POINTER 3
#define ALIB_ERROR_ARRAY_SIZE 4
#define ALIB_ERROR_OUT_OF_BOUND 5

#define ALIB_Error(c) 

typedef struct {
	double* frequency;
	double min_range;
	double max_range;
	double d_range;
	unsigned int length;
	unsigned int samples;
} alib_histogram;

typedef struct {
	double *y_data;
	double dx;  //usually the min_x/max_x are not so interesting anyway
	double xshift;
	unsigned int length;
	unsigned int samples; //for calculating averages
} alib_function;

typedef alib_function alib_tensor;


alib_histogram* ALIB_HistogramMalloc(double min_r, double max_r, unsigned int length);

int ALIB_HistogramFree(alib_histogram* p_ahisto);

alib_function* ALIB_FunctionMalloc(double dx, double xshift, unsigned int length);

int ALIB_FunctionFree(alib_function* p_afxn);


//int ALIB_FunctionDivid(alib_function* p_afxn,alib_function* p_bfxn); 

int  ALIB_Correlation( unsigned short mode, alib_function* p_afxn, 
		double* sample_data,unsigned int dim, unsigned int length,
		double (*corr_fxn)(double* ,double*,unsigned int),FILE* f_out); //in diffusion this would calculate <r(0).r(t)>

int ALIB_Histogram(unsigned short mode, alib_histogram* p_ahisto, 
		double* sample_data,unsigned int dim, unsigned int length,double (*_fxn)(double*,unsigned int),FILE* f_out);

int ALIB_AverageProfileFunction(unsigned short mode, alib_function *p_afxn, alib_function* p_afreq,
		double* xshift, double* sample_data,unsigned int length, FILE *f_out);

#define ALIB_TensorMalloc(rank) ALIB_FunctionMalloc((double) rank,0.0, (unsigned int) rank*rank)
//instead of rank*rank, it's probably better to use rank*(rank+1)/2, uses less
//memory
#define ALIB_TensorFree(p_tensor) ALIB_FunctionFree(p_tensor)
int ALIB_InertiaTensor(unsigned short mode, alib_tensor* p_tensor,
		double* mass_and_pos, unsigned int dim_plus_1, unsigned int length,
		FILE* f_out);


#endif
