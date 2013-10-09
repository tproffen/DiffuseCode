#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cuComplex.h"

#define I2PI (1<<16)
#define MASK (I2PI-1)

__global__ void init_array_kernel(cuFloatComplex*, int);

__global__ void four_strucf_cuda_kernel(cuFloatComplex*, cuFloatComplex*,
					int, int, int,
					int, int, int);

extern "C"{
  void four_strucf_cuda_(cuFloatComplex *tcsf, cuFloatComplex *cex, float *xat, int *nxat, int *num, float *xm, float *win, float *vin, float *uin, int *cr_natoms)
  {
    int nnum = num[0]*num[1]*num[2];
    double xarg0, xincu, xincv;
    int iarg0, iincu, iincv;
    int blockDim = 64;
    int gridDim = (nnum + blockDim - 1) / blockDim;
    
    //Allocate space for complex exponent and copy to cex to device.
    cuFloatComplex* d_cex;
    cudaMalloc((void**) &d_cex, I2PI * sizeof(cuFloatComplex));
    cudaMemcpy(d_cex, cex, I2PI * sizeof(cuFloatComplex), cudaMemcpyHostToDevice);
    
    //Allocate space for csf table and initalise to zero.
    cuFloatComplex* d_tcsf;
    cudaMalloc((void**) &d_tcsf, nnum * sizeof(cuFloatComplex));
    init_array_kernel<<<gridDim, blockDim>>>(d_tcsf, nnum);
    
    //Loop over all of the atoms we are handling now ... 
    for(int l=0; l< nxat[0]; l++){
      // Get initial argument to the exponent and increments along the two axies 'u' and 'v'.
      xarg0 = xm[0]  * xat[l] + xm[1]  * xat[l+cr_natoms[0]] + xm[2]  * xat[l+cr_natoms[0]*2];
      xincu = uin[0] * xat[l] + uin[1] * xat[l+cr_natoms[0]] + uin[2] * xat[l+cr_natoms[0]*2];
      xincv = vin[0] * xat[l] + vin[1] * xat[l+cr_natoms[0]] + vin[2] * xat[l+cr_natoms[0]*2];
      //Convert to high precision integers (64*i2pi=2^20) ...
      iarg0 = (int) rint( 64 * I2PI * (xarg0 - (int) xarg0 + 1.) );
      iincu = (int) rint( 64 * I2PI * (xincu - (int) xincu + 1.) );
      iincv = (int) rint( 64 * I2PI * (xincv - (int) xincv + 1.) );
      
      four_strucf_cuda_kernel<<<gridDim, blockDim>>>
	(d_cex, d_tcsf,
	 num[0],num[1],nnum,
	 iarg0,iincu,iincv);
    }
    
    cudaMemcpy(tcsf, d_tcsf, nnum*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost);
    
    cudaFree(d_tcsf);
    cudaFree(d_cex);
    
  }
}

__global__ void four_strucf_cuda_kernel(cuFloatComplex* d_cex, cuFloatComplex* d_tcsf,
					int num1, int num2, int nnum,
					int iarg0, int iincu, int iincv)
{
  int i, j, iadd, iarg;
  unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x;
  //Loop over all image pixels 'idx'. 'iadd' is the address of the
  //argument to the complex exponent (in the table 'cex()'). The '>>6'
  //operation divides out the 64 and the '&MASK' is used so that the
  //argument to the complex exponent is inside our table which has
  //range 0=>2pi.
  if(idx<nnum)
    {
      i = idx / num1;
      j = idx % num1;
      iarg = iarg0 + i * iincu + j * iincv;
      iadd = iarg >> 6;
      iadd = iadd & MASK;
      d_tcsf[idx] = cuCaddf(d_tcsf[idx],d_cex[iadd]);
    };
  __syncthreads();
}

__global__ void init_array_kernel(cuFloatComplex* array, int n)
{
  unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x;
  if(idx<n)
    array[idx] = make_cuFloatComplex(0.0,0.0);
  __syncthreads();
}
