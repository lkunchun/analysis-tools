#include "analysis_lib.h"
#define MAX_DATA_POINTS 24000
#define SEED 123456
#define DIM_PLUS_ONE 4
#define CYLINDER_RADIUS 0.1
#define CYLINDER_FREQUENCY 1000
#define CYLINDER_DHEIGHT 0.001

double rtn_double(double* p_double,unsigned int d) {
	return *p_double;
}

int main(void) {

	double test_data[MAX_DATA_POINTS];
	int i_count;
	alib_histogram* p_hist;
	alib_tensor* p_tensor;

	srand48(SEED);
	for(i_count=0;i_count<MAX_DATA_POINTS;i_count++) {
		test_data[i_count]=drand48();
	}


	//testing histogram routine
	
	p_hist = ALIB_HistogramMalloc(0.0,1.0,100);
	ALIB_Histogram(ALIB_SAMPLING_MODE,p_hist,test_data,1,MAX_DATA_POINTS, rtn_double,NULL);
	ALIB_Histogram(ALIB_PRINT_MODE,p_hist,NULL,1,MAX_DATA_POINTS, rtn_double,stderr);

	//////////////////////////////////////////////////////////////////////////
	//
	//testing inertial tensor routine
	
	p_tensor = ALIB_TensorMalloc(3);
	for(i_count=0;i_count<MAX_DATA_POINTS;i_count++) {
		if(i_count%DIM_PLUS_ONE==0) test_data[i_count]=1.0/MAX_DATA_POINTS*DIM_PLUS_ONE;
		else test_data[i_count]-=0.5;
	}
	fprintf(stderr,"#inertia tensor for random mass points\n");
	ALIB_InertiaTensor(ALIB_SAMPLING_MODE,p_tensor,test_data,DIM_PLUS_ONE,MAX_DATA_POINTS,NULL);
	ALIB_InertiaTensor(ALIB_PRINT_MODE,p_tensor,NULL,DIM_PLUS_ONE,MAX_DATA_POINTS,stderr);
	//for isotropic distribution of masses in three D the off diagonal
	//entries should be 0, while the main diagonal is (dim-1)*1/6
	//
	//next try cylindrical distribution with symmetry axes in z direction
	//assume DIM_PLUS_ONE == 4
	if(DIM_PLUS_ONE==4) {
		ALIB_InertiaTensor(ALIB_INITIALIZATION_MODE,p_tensor,NULL,DIM_PLUS_ONE,MAX_DATA_POINTS,NULL);
		for(i_count=0;i_count<MAX_DATA_POINTS;i_count+=DIM_PLUS_ONE) {
			test_data[i_count]=1.0/MAX_DATA_POINTS*DIM_PLUS_ONE;
			test_data[i_count+1] = CYLINDER_RADIUS*cos(M_PI*2.0/CYLINDER_FREQUENCY*i_count/DIM_PLUS_ONE);
			test_data[i_count+2] = CYLINDER_RADIUS*sin(M_PI*2.0/CYLINDER_FREQUENCY*i_count/DIM_PLUS_ONE);
			test_data[i_count+3] = (CYLINDER_DHEIGHT)*(((float) i_count)/DIM_PLUS_ONE);
		}
		fprintf(stderr,"#inertia tensor for cylindrical mass points\n");
		ALIB_InertiaTensor(ALIB_SAMPLING_MODE,p_tensor,test_data,DIM_PLUS_ONE,MAX_DATA_POINTS,NULL);
		ALIB_InertiaTensor(ALIB_PRINT_MODE,p_tensor,NULL,DIM_PLUS_ONE,MAX_DATA_POINTS,stderr);
	}

	return ALIB_SUCCESS;
}
