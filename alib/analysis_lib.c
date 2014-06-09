#include "analysis_lib.h"


alib_histogram* ALIB_HistogramMalloc(double min_r, double max_r, unsigned int length) {
	alib_histogram* p_rtn=NULL;
	int i_count;
	p_rtn=malloc(sizeof(alib_histogram));
	if(p_rtn==NULL) {
		ALIB_Error(ALIB_ERROR_NULL_POINTER);
		return NULL;
	}
	p_rtn->min_range=(min_r<max_r)?min_r:max_r;
	p_rtn->max_range=(min_r<max_r)?max_r:min_r;
	p_rtn->length=length;
	p_rtn->samples = 0;
	p_rtn->d_range=fabs(max_r-min_r)/length;
	p_rtn->frequency=malloc(length*sizeof(double));
	if(p_rtn->frequency==NULL) {
		ALIB_Error(ALIB_ERROR_NULL_POINTER);
		return NULL;
	}
	for(i_count=0;i_count<length;i_count++) {
		(p_rtn->frequency)[i_count]=0.0;
	}

	return p_rtn;
}

int ALIB_HistogramFree(alib_histogram* p_ahisto) {
	if(p_ahisto!=NULL) {
		free(p_ahisto->frequency);
		free(p_ahisto);
	} else {
		ALIB_Error(ALIB_ERROR_NULL_POINTER);
		return ALIB_ERROR;
	}
	return ALIB_SUCCESS;
}

alib_function* ALIB_FunctionMalloc(double dx, double xshift, unsigned int length) {
	alib_function* p_rtn=NULL;
	int i_count;
	p_rtn=malloc(sizeof(alib_function));
	if(p_rtn==NULL) {
		ALIB_Error(ALIB_ERROR_NULL_POINTER);
		return NULL;
	}
	p_rtn->dx = dx;
	p_rtn->length=length;
	p_rtn->xshift=xshift;
	p_rtn->samples = 0;
	p_rtn->y_data=malloc(length*sizeof(double));
	if(p_rtn->y_data==NULL) {
		ALIB_Error(ALIB_ERROR_NULL_POINTER);
		return NULL;
	}
	for(i_count=0;i_count<length;i_count++) {
		(p_rtn->y_data)[i_count]=0.0;
	}
	return p_rtn;
}


int ALIB_FunctionFree(alib_function* p_afxn) {
	if(p_afxn!=NULL) {
		free(p_afxn->y_data);
		free(p_afxn);
	} else {
		ALIB_Error(ALIB_ERROR_NULL_POINTER);
		return ALIB_ERROR;
	}
	return ALIB_SUCCESS;
}


int ALIB_AverageProfileFunction(unsigned short mode, alib_function *p_afxn, alib_function *p_afreq,
		double *xshift,
		double* sample_data,unsigned int data_length, 
		FILE *f_out) {

	int i_count,k_count;
	double result;
	double* p_double;

	if(p_afxn==NULL) {
		ALIB_Error(ALIB_ERROR_NULL_POINTER);
		return ALIB_ERROR;
	}

	switch(mode) {

		case ALIB_INITIALIZATION_MODE:  //user should just kill p_ahisto and restart another handler
			p_afxn->samples=0;
			p_afreq->samples=0;
			for(i_count=0;i_count< (p_afxn->length);i_count++) {
				(p_afxn->y_data)[i_count]=0.0;
			}
			for(i_count=0;i_count< (p_afreq->length);i_count++) {
				(p_afreq->y_data)[i_count]=0.0;
			}

			if(p_afreq->length<p_afxn->length) {
				fprintf(stderr,"ALIB error: incompatible size\n");
				return  ALIB_ERROR;
			}
			break;

			//xshift won't get reset

		case ALIB_SAMPLING_MODE:
			for(k_count=0;k_count<data_length;k_count++) {
				i_count=(int) ((xshift[k_count]-p_afxn->xshift)/p_afxn->dx);
				if(i_count<0 || i_count>=p_afxn->length) {
					fprintf(stderr,"ALIB error: negative index\n");
					ALIB_Error(ALIB_ERROR_OUT_OF_BOUND);
					return ALIB_ERROR;
				} else {
					result=sample_data[k_count];
					(p_afxn->y_data)[i_count]+=result;
					
					(p_afreq->y_data)[i_count]+=1.0;
				}
			}
			p_afxn->samples++;
			break;

		case ALIB_PRINT_MODE:

			if(p_afxn->samples==0) return ALIB_ERROR;
			//result = (xshift-p_afxn->xshift);
			result=0.0;
			fprintf(f_out,"#alib profile function\n");
			for(i_count=0;i_count<p_afxn->length;i_count++) {
				
				//if(data_length==0) {
				//	fprintf(f_out,"%g %g\n",result,(p_afxn->y_data)[i_count]);
				//} else {
				//	fprintf(f_out,"%g %g\n",result,(p_afxn->y_data)[i_count]/data_length);
				//}
				//

				fprintf(f_out,"%g %g\n",result,(p_afxn->y_data)[i_count]/(p_afreq->y_data)[i_count]);
				result+=p_afxn->dx;

			}
			fprintf(f_out,"#end profile function\n");
			break;

		default:
			fprintf(stderr,"unrecognized ALIB token\n");
			exit(0); //should really learn how to throw signal
			break;

	}

	return ALIB_SUCCESS;
}

int ALIB_Histogram(unsigned short mode, alib_histogram* p_ahisto, 
		double* sample_data,unsigned int dim, unsigned int length,double (*_fxn)(double*,unsigned int), FILE* f_out) {

	int i_count;
	double result;
	double* p_double;

	if(p_ahisto==NULL) {
		ALIB_Error(ALIB_ERROR_NULL_POINTER);
		return ALIB_ERROR;
	}

	switch (mode) {
		case ALIB_INITIALIZATION_MODE:  //user should just kill p_ahisto and restart another handler

			p_ahisto->samples = 0;

			for(i_count=0;i_count<p_ahisto->length;i_count++) {
				(p_ahisto->frequency)[i_count]=0.0;
			}

			break;



		case ALIB_SAMPLING_MODE:


			p_ahisto->samples++;


			for(i_count=0;i_count<length;i_count+=dim) {

				p_double=(sample_data)+i_count;
				result=_fxn(p_double,dim);
				if(result>(p_ahisto->min_range) && result < (p_ahisto->max_range)) {
					result-=p_ahisto->min_range;
					result/=p_ahisto->d_range;
				}
				(p_ahisto->frequency)[(int) result]++; //left boundary

			}

			break;

		case ALIB_PRINT_MODE:
			if(p_ahisto->samples==0) return ALIB_ERROR;
			result = p_ahisto->min_range;
			fprintf(f_out,"#alib histogram\n");
			for(i_count=0;i_count<p_ahisto->length;i_count++) {
				fprintf(f_out,"%g %g\n",result,(p_ahisto->frequency)[i_count]/(p_ahisto->samples));
				result+=p_ahisto->d_range;
			}
			fprintf(f_out,"#end histogram\n");
			break;
		default:
			fprintf(stderr,"unrecognized ALIB token\n");
			exit(0); //should really learn how to throw signal
			break;
	}
	return ALIB_SUCCESS;

}

int  ALIB_Correlation( unsigned short mode, alib_function* p_afxn, 
		double* sample_data,unsigned int dim, unsigned int length,
		double (*corr_fxn)(double* ,double*,unsigned int ),FILE* f_out) {//in diffusion this would calculate <r(0).r(t)>

	int i_count, j_count;
	double result;
	double *p_double1, *p_double2;

	if(p_afxn==NULL) {
		ALIB_Error(ALIB_ERROR_NULL_POINTER);
		return ALIB_ERROR;
	}
	
	switch (mode) {

		case ALIB_INITIALIZATION_MODE:  

			p_afxn->samples = 0;

			for(i_count=0;i_count<p_afxn->length;i_count++) {
				(p_afxn->y_data)[i_count]=0.0;
			}

			break;

		case ALIB_SAMPLING_MODE:

			for(i_count=0;i_count<(length-p_afxn->length*dim);i_count+=dim) {
				p_afxn->samples++;
				p_double1=sample_data+i_count;
				for(j_count=i_count;j_count<(i_count+p_afxn->length*dim);j_count+dim) {
					p_double2=sample_data+j_count;
					(p_afxn->y_data)[j_count/dim]+=corr_fxn(p_double1,p_double2,dim);
				}
			}

			break;

		case ALIB_PRINT_MODE:
			if(p_afxn->samples==0) return ALIB_ERROR;
			result = p_afxn->xshift;
			fprintf(f_out,"#alib correlation function\n");
			for(i_count=0;i_count<p_afxn->length;i_count++) {
				fprintf(f_out,"%g %g\n",result,(p_afxn->y_data)[i_count]/(p_afxn->samples));
				result+=p_afxn->dx;
			}
			fprintf(f_out,"#end correlation function\n");
			break;
		default:
			fprintf(stderr,"unrecognized ALIB token\n");
			exit(0); //should really learn how to throw signal
			break;
	}

	return ALIB_SUCCESS;

}

int ALIB_InertiaTensor(unsigned short mode, alib_tensor* p_tensor,
		double* mass_and_pos, unsigned int dim_plus_1, unsigned int length,
		FILE* f_out) {

	int i_count,j_count,k_count, dim;
	double *p_double; 
	double r2;

	if(p_tensor==NULL) {
		ALIB_Error(ALIB_ERROR_NULL_POINTER);
		return ALIB_ERROR;
	}
	if(p_tensor->length != (dim_plus_1-1)*(dim_plus_1-1)) {
		ALIB_Error(ALIB_ERROR_ARRAY_SIZE);
		return ALIB_ERROR;
	}

	dim = dim_plus_1-1;
	
	switch (mode) {

		case ALIB_INITIALIZATION_MODE:  

			p_tensor->samples = 0;

			for(i_count=0;i_count<p_tensor->length;i_count++) {
				(p_tensor->y_data)[i_count]=0.0;
			}

			break;

		case ALIB_SAMPLING_MODE:
			//see page 408 Marion, "Classical Dynamics"
			
			p_tensor->samples++;
			for(i_count=0;i_count<length;i_count+=dim_plus_1) {

				p_double=mass_and_pos+i_count;
				r2=0.0;

				for(j_count=1;j_count<dim_plus_1;j_count++) {
					r2+=p_double[j_count]*p_double[j_count];
				}

				r2*=(*p_double);
				
				for(j_count=0;j_count<dim;j_count++) {
					for(k_count=0;k_count<dim;k_count++) {
						(p_tensor->y_data)[j_count+k_count*dim]
							-= (*p_double)*(*(p_double+j_count+1))*(*(p_double+k_count+1));
						if(j_count==k_count) 
							(p_tensor->y_data)[j_count+k_count*dim]+=r2;
					}
				}
			}

			break;

		case ALIB_PRINT_MODE:
			if(p_tensor->samples==0) return ALIB_ERROR;
			fprintf(f_out,"#alib inertial tensor\n");
			for(i_count=0;i_count<sqrt(p_tensor->length);i_count++) {
				for(j_count=0;j_count<sqrt(p_tensor->length);j_count++) {
					fprintf(f_out," %f",
						(p_tensor->y_data)[i_count+j_count*(dim_plus_1-1)]/(p_tensor->samples));
				}
				fprintf(f_out,"\n");
			}
			fprintf(f_out,"#end inertia tensor\n");
			break;
		default:
			fprintf(stderr,"unrecognized ALIB token\n");
			exit(0); //should really learn how to throw signal
			break;
	}

	return ALIB_SUCCESS;
}


