#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cuComplex.h"

#define I2PI (1<<16)
#define MASK (I2PI-1)

__global__ void init_array_kernel(cuFloatComplex*, int);

__global__ void four_strucf_cuda_kernel(cuFloatComplex*, cuFloatComplex*,
					int, int, int, int, int,
					int*, int*, int*, int*);

__global__ void iin_xat(float*, int*, int*, int*, int*,
			float, float, float,
			float, float, float,
			float, float, float,
			float, float, float,
			int, int);

extern "C"{
  void four_strucf_cuda_(cuFloatComplex *tcsf, cuFloatComplex *cex,
			 float *xat, int *nxat, int *num,
			 float *xm, float *win, float *vin, float *uin, int *cr_natoms)
  {
    int nnum = num[0]*num[1]*num[2];
    int blockDim = 64;
    int gridDim = (nnum + blockDim - 1) / blockDim;
    int gridDim2 = (cr_natoms[0]+blockDim - 1 ) / blockDim;
    
    //Allocate space for complex exponent and copy to cex to device.
    cuFloatComplex* d_cex;
    cudaMalloc((void**) &d_cex, I2PI * sizeof(cuFloatComplex));
    cudaMemcpy(d_cex,      cex, I2PI * sizeof(cuFloatComplex), cudaMemcpyHostToDevice);
    
    //Allocate space for csf table and initalise to zero.
    cuFloatComplex* d_tcsf;
    cudaMalloc((void**) &d_tcsf, nnum * sizeof(cuFloatComplex));
    init_array_kernel<<<gridDim, blockDim>>>(d_tcsf, nnum);
    
    //Allocate space for iarg0, iincu, iincv, iincw
    int* d_iarg0;
    int* d_iincu;
    int* d_iincv;
    int* d_iincw;
    cudaMalloc((void**) &d_iarg0, nxat[0] * sizeof(int));
    cudaMalloc((void**) &d_iincu, nxat[0] * sizeof(int));
    cudaMalloc((void**) &d_iincv, nxat[0] * sizeof(int));
    cudaMalloc((void**) &d_iincw, nxat[0] * sizeof(int));
    
    //Allocate space for list of atoms (xat) and copy to device
    float* d_xat;
    cudaMalloc((void**) &d_xat, cr_natoms[0] * 3 * sizeof(float));
    cudaMemcpy(d_xat,      xat, cr_natoms[0] * 3 * sizeof(float), cudaMemcpyHostToDevice);
    
    iin_xat<<<gridDim2, blockDim>>>
      (d_xat, d_iarg0, d_iincu, d_iincv, d_iincw,
       xm[0], xm[1], xm[2],
       uin[0], uin[1], uin[2],
       vin[0], vin[1], vin[2],
       win[0], win[1], win[2],
       cr_natoms[0], nxat[0]);
    
    four_strucf_cuda_kernel<<<gridDim, blockDim>>>
      (d_cex, d_tcsf,
       num[0],num[1],num[2],nnum,nxat[0],
       d_iarg0, d_iincu, d_iincv, d_iincw);
    
    cudaMemcpy(tcsf, d_tcsf, nnum*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost);
    
    cudaFree(d_tcsf);
    cudaFree(d_cex);
    cudaFree(d_xat);
    cudaFree(d_iarg0);
    cudaFree(d_iincu);
    cudaFree(d_iincv);
    cudaFree(d_iincw);
    
  }
}

__global__ void iin_xat(float* d_xat, int* d_iarg0, int* d_iincu, int* d_iincv, int* d_iincw,
			float xm1,  float xm2,  float xm3,
			float uin1, float uin2, float uin3,
			float vin1, float vin2, float vin3,
			float win1, float win2, float win3,
			int cr_natoms, int nxat)
{
  unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x;
  double xarg0, xincu, xincv, xincw;
  
  if(idx<nxat)
    {
      xarg0 = xm1  * d_xat[idx] + xm2  * d_xat[idx+cr_natoms] + xm3  * d_xat[idx+cr_natoms*2];
      xincu = uin1 * d_xat[idx] + uin2 * d_xat[idx+cr_natoms] + uin3 * d_xat[idx+cr_natoms*2];
      xincv = vin1 * d_xat[idx] + vin2 * d_xat[idx+cr_natoms] + vin3 * d_xat[idx+cr_natoms*2];
      xincw = win1 * d_xat[idx] + win2 * d_xat[idx+cr_natoms] + win3 * d_xat[idx+cr_natoms*2];
      d_iarg0[idx] = (int)rintf(64 * I2PI * (xarg0 - (int)xarg0 + 1.));
      d_iincu[idx] = (int)rintf(64 * I2PI * (xincu - (int)xincu + 1.));
      d_iincv[idx] = (int)rintf(64 * I2PI * (xincv - (int)xincv + 1.));
      d_iincw[idx] = (int)rintf(64 * I2PI * (xincw - (int)xincw + 1.));
    }
}

__global__ void four_strucf_cuda_kernel(cuFloatComplex* d_cex, cuFloatComplex* d_tcsf,
					int num1,   int   num2, int num3, int nnum, int nxat,
					int* d_iarg0, int* d_iincu, int* d_iincv, int* d_iincw)
{
  int i, j, k, n, iadd, iarg;
  unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x;
  //Loop over all image pixels 'idx'. 'iadd' is the address of the
  //argument to the complex exponent (in the table 'cex()'). The '>>6'
  //operation divides out the 64 and the '&MASK' is used so that the
  //argument to the complex exponent is inside our table which has
  //range 0=>2pi.
  
  if(idx<nnum)
    {
      i = idx / (num2 * num3);
      j = (idx / num3) % num2;
      k = idx % num3;
      for(n=0; n < nxat; n++)
	{
	  iarg = d_iarg0[n] + i * d_iincu[n] + j * d_iincv[n] + k * d_iincw[n];
	  iadd = iarg >> 6;
	  iadd = iadd & MASK;
	  d_tcsf[idx] = cuCaddf(d_tcsf[idx],d_cex[iadd]);
	}
    }
  __syncthreads();
}

__global__ void init_array_kernel(cuFloatComplex* array, int n)
{
  unsigned int idx = blockDim.x*blockIdx.x + threadIdx.x;
  if(idx<n)
    array[idx] = make_cuFloatComplex(0.0,0.0);
  __syncthreads();
}
