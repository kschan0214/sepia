/*  mexUnwrap.c
    --------------------
    written by Tian Liu

    --------------------
    interfaced with matlab mex by Dong Zhou
    zhou.dong@gmail.com
    ---------------------
    created:     7.10.13
    last update: 7.11.13

*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mex.h"

#define PI 3.1415926535897931

struct MRI_point {
    int pos[3];
	int ind;
	float mag;
	float phase;
	int planned;
	int pre;
};
typedef struct MRI_point point;

static float round (float num){
    return floor(num+0.5);
}
static void swap(int *a, int *b)
{
int c;
	c = *a;
	*a = *b;
	*b = c;
}

static void get_neighbors(int *nei,int *pos,int *dim)
{
    
    int x,y,z;

	x = pos[0];
	y = pos[1];
	z = pos[2];

	if(x>0)
        nei[0]=x-1+y*dim[0]+z*dim[0]*dim[1];
    else
        nei[0]=-1;
    
    if (x<(dim[0]-1))
        nei[1]=x+1+y*dim[0]+z*dim[0]*dim[1];
    else
        nei[1]=-1;

    if (y>0)
        nei[2]=x+(y-1)*dim[0]+z*dim[0]*dim[1];
    else
        nei[2]=-1;
    
    if (y<(dim[1]-1))
        nei[3]=x+(1+y)*dim[0]+z*dim[0]*dim[1];
    else
        nei[3]=-1;

    if (z>0)
        nei[4]=x+y*dim[0]+(z-1)*dim[0]*dim[1];
    else
        nei[4]=-1;
    
    if (z<(dim[2]-1))
        nei[5]=x+y*dim[0]+(z+1)*dim[0]*dim[1];
    else
        nei[5]=-1;

}

static point* assemble_data( float *mag, float *phase, int *dim, int *max_ind)
{
int ind, numel;

point * ret_image;

	numel = dim[0]*dim[1]*dim[2];
	ret_image = (point *) calloc(numel, sizeof(point));

	for (ind=0; ind< numel; ind++)
	{
    	ret_image[ind].pos[2]=(int)((ind)/(dim[1]*dim[0]));
    	ret_image[ind].pos[1]=(int)((ind - ret_image[ind].pos[2] *dim[1]*dim[0])/dim[0]);
    	ret_image[ind].pos[0]=(int)(ind - ret_image[ind].pos[1]*dim[0] - ret_image[ind].pos[2]*dim[1]*dim[0]);
		ret_image[ind].ind = ind;
		ret_image[ind].mag = mag[ind];
		ret_image[ind].phase = phase[ind];
		ret_image[ind].planned = 0;
		ret_image[ind].pre = 0;
		if (mag[ind] > mag[*max_ind]) *max_ind = ind;
	}	
	
	return ret_image;

}

static int * heap_create(int seed, int * dim)
{
int numel;
int *heap;
	numel = dim[0]*dim[1]*dim[2];
	heap = (int *) calloc(numel, sizeof(int));
	heap[0] = seed;
	return heap;
}

static void heap_add(point *ori_image, int * heap, int index, int * length)
{
int ind, parent, goon;
	
	heap[*length] = index;
	goon = 1;
	
	ind = *length;
	while ((ind>0) && (goon))
	{
		parent = (ind-1)/2;
	
		if (ori_image[heap[parent]].mag < ori_image[heap[ind]].mag) 
		{
			swap( &heap[parent], &heap[ind]);
			ind = parent;
		}
		else
		{
			goon = 0;
		}
	}
	*length = *length + 1;
}

static int heap_pop( point * ori_image, int * heap, int * length)
{
int left, right, ind, goon, child;
int ret_val;	
	*length = *length - 1;
	ret_val = heap[0];
	heap[0] = heap[*length];
	heap[*length] = 0;
	
	ind = 0;
	goon = 1;
	while ((ind*2+1<= *length - 1) && (goon))
	{
		left = ind * 2 +1;
		right = ind*2 +2;
		if ( ori_image[heap[left]].mag > ori_image[heap[right]].mag )
			child = left;
		else
			child = right;
		if (child > (*length - 1)) child = left;
		
		if (ori_image[heap[ind]].mag < ori_image[heap[child]].mag)
		{
			swap( &heap[ind], &heap[child]);
			ind = child;
		}
		else
		{
			goon = 0;
		}
	}
	return ret_val;
}

void unwrap(point *ori_image, int seed, int *dim, float noise)
{
int nei[6];
int curr_ind, counter;
int i;
int * heap;
int heap_length = 0;

	counter = 0;
	curr_ind = seed;

	/* beginning of the searching algorithm	*/
    
	heap = heap_create(seed, dim);
	heap_length++;
	ori_image[heap[0]].pre = -1;
	ori_image[heap[0]].planned = 1;
		
	do
	{

		counter ++; 
/*		
		if (ori_image[heap[0]].phase < 10)
			ori_image[heap[0]].phase = 10;
		else 
			ori_image[heap[0]].phase = ori_image[heap[0]].phase + 10;
*/		
		
		if (ori_image[heap[0]].pre >=0) 
			ori_image[heap[0]].phase = ori_image[heap[0]].phase - round( ( ori_image[heap[0]].phase - ori_image[ ori_image[heap[0]].pre ].phase )/(2*PI) )   *2*PI;

		get_neighbors(nei, ori_image[heap[0]].pos, dim);

		curr_ind = heap_pop( ori_image, heap, &heap_length);

        
		for (i = 0; i<6; i++)
        {
            if ( (nei[i]>=0) && (0 == ori_image[nei[i]].planned ) && (ori_image[nei[i]].mag > noise) )
			{
				ori_image[nei[i]].pre = curr_ind;
				heap_add(ori_image, heap, nei[i], &heap_length);
				
				ori_image[nei[i]].planned = 1;
			}
        }

		
	} while (heap_length>0 );

	printf("%d voxels are unwrapped\n", counter);
	free(heap); heap = NULL;
    
}

/* interface */
void mexFunction( int nlhs, mxArray      *plhs[],
		          int nrhs, const mxArray *prhs[]) {
    int dim[3], i;
    mxArray *out_m;
    int numel;
    double *c;
    float noise_level = 0;
    float noise_ratio = 5e-3;
    point * MRimage;
    float *phase, *mag;
	int max_ind=0;	/* the index of the element with largest intensity */
    FILE* fid_in;
    double *p;

    if (nrhs != 5) {
    mexErrMsgTxt("MEX requires 5 input arguments: magnitude_map phase_map dimension1 dimension2 dimension3"); }
    
    dim[0] = (int) (*mxGetPr(prhs[2]));
    dim[1] = (int) (*mxGetPr(prhs[3]));
    dim[2] = (int) (*mxGetPr(prhs[4]));
    numel = dim[0] * dim[1] * dim[2];

    mag = (float *) calloc( numel, sizeof(float));
    phase = (float *) calloc( numel, sizeof(float));

    /* read data */ 
    p = mxGetPr(prhs[0]);
    for (i=0;i<numel;++i){
        mag[i] = p[i];
    }
    p = mxGetPr(prhs[1]);
    for (i=0;i<numel;++i){
        phase[i] = p[i];
    }

    /* calculate */
    MRimage = assemble_data(mag,phase,dim,&max_ind);
	noise_level = MRimage[max_ind].mag * noise_ratio;

	unwrap(MRimage, dim[2]/2*dim[1]*dim[0] + dim[1]/2*dim[0] + dim[0]/2,  dim, noise_level);	/* choose the middle of the FOV as the seeding point	*/

    out_m = plhs[0] = mxCreateDoubleMatrix(numel,1,mxREAL);

    c = mxGetPr(out_m);
    for (i=0;i<numel;++i)
        c[i] = MRimage[i].phase;

	free(MRimage); 
	free(mag); 
	free(phase); 
} 





