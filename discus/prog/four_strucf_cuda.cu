#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cuComplex.h"

#define I2PI (1<<16)
#define MASK (I2PI-1)

__global__ void initarraygpu(cuFloatComplex[], int);

__global__ void computestrucf(cuFloatComplex*, cuFloatComplex*,
			      int, int,
			      int, int, int);

extern "C"{
  void cudastrucf_(cuFloatComplex *csf, cuFloatComplex *cex, float *xat, int *nxat, int *num, float *xm, float *win, float *vin, float *uin, int *cr_natoms)
  {
    int nnum = num[0]*num[1]*num[2];
    
    int threadsPerBlock = 64;
    int threadsPerGrid = (nnum + threadsPerBlock - 1) / threadsPerBlock;
    
    cuFloatComplex* d_tcsf;
    cudaMalloc((void**) &d_tcsf, nnum * sizeof(cuFloatComplex));
    
    cuFloatComplex* d_cex;
    cudaMalloc((void**) &d_cex, I2PI * sizeof(cuFloatComplex));
    
    cudaMemcpy(d_cex, cex, I2PI * sizeof(cuFloatComplex), cudaMemcpyHostToDevice);
    
    initarraygpu<<<threadsPerGrid, threadsPerBlock>>>(d_tcsf, nnum);
    
    printf("Starting CUDA!!!!\n");
    
    float xarg0, xincu, xincv;//, xincw;
    int iarg0, iincu, iincv;//, iincw;
    
    for(int l=0; l< nxat[0]; l++){
      xarg0 = xm[0] * xat[l] + xm[1] * xat[l+cr_natoms[0]+1] + xm[2] * xat[l+cr_natoms[0]+2];
      xincu = uin[0] * xat[l] + uin[1] * xat[l+cr_natoms[0]+1] + uin[2] * xat[l+cr_natoms[0]+2];
      xincv = vin[0] * xat[l] + vin[1] * xat[l+cr_natoms[0]+1] + vin[2] * xat[l+cr_natoms[0]+2];
      //xincw = win1 * xat1 + win2 * xat2 + win3 * xat3;
      iarg0 = (int)rintf(64 * I2PI * (xarg0 - (int)xarg0 + 1.));
      iincu = (int)rintf(64 * I2PI * (xincu - (int)xincu + 1.));
      iincv = (int)rintf(64 * I2PI * (xincv - (int)xincv + 1.));
      //iincw = (int)rintf(64 * I2PI * (xincw - (int)xincw + 1.));
      
      computestrucf<<<threadsPerGrid, threadsPerBlock>>>
	(d_cex, d_tcsf,
	 num[0],num[1],
	 iarg0,iincu,iincv);
    }
    
    cudaMemcpy(csf, d_tcsf, nnum*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost);
    cudaMemcpy(csf, d_tcsf, nnum*sizeof(cuFloatComplex), cudaMemcpyDeviceToHost);
        
    cudaFree(d_tcsf);
    cudaFree(d_cex);
    
  }
}

__global__ void computestrucf(cuFloatComplex* cex, cuFloatComplex* tcsf,
			      int num1, int num2,
			      int iarg0, int iincu, int iincv)
{
  int i, j, iadd, id, iarg;
  
  id = threadIdx.x + blockDim.x * blockIdx.x;
  if(id<num1*num2)
    {
      i = id / num1;
      j = id % num1;
      iarg = iarg0 + i * iincu + j * iincv;
      iadd = iarg >> 6;
      iadd = iadd & MASK;
      tcsf[id] = cuCaddf(tcsf[id],cex[iadd]);
    };
  __syncthreads();
}


__global__ void initarraygpu(cuFloatComplex* array1, int nelements)
{
  int id = threadIdx.x + blockDim.x * blockIdx.x;
  if(id<nelements)
    {
      array1[id] = make_cuFloatComplex(0.0,0.0);
    };
  __syncthreads();
}

