#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define I2PI (1<<16)
#define MASK (I2PI-1)

__global__ void initarraygpu(float[], float[], int);

__global__ void computestrucf(float, float, float, float[], float[], float[], float[],
			      int, int, float, float, float,
			      float, float, float, float, float, float, 
			      float, float, float);

extern "C" int __config_mod_MOD_nmax;
extern "C" float __diffuse_mod_MOD_xm[3], __diffuse_mod_MOD_win[3], __diffuse_mod_MOD_vin[3], __diffuse_mod_MOD_uin[3];
extern "C" int  __diffuse_mod_MOD_num[3], __diffuse_mod_MOD_nxat;
extern "C" int __crystal_mod_MOD_cr_natoms;

extern "C"{
  void cudastrucf_(float *csf_r, float *csf_i, float *cex_r, float *cex_i, float *xat)
  {
    int nnum = __diffuse_mod_MOD_num[0]*__diffuse_mod_MOD_num[1]*__diffuse_mod_MOD_num[2];
    
    int threadsPerBlock = 64;
    int threadsPerGrid = (nnum + threadsPerBlock - 1) / threadsPerBlock;
    
    float* d_rtcsf;
    cudaMalloc((void**) &d_rtcsf, nnum * sizeof(float));
    float* d_itcsf;
    cudaMalloc((void**) &d_itcsf, nnum * sizeof(float));
    
    float* d_rexp;
    cudaMalloc((void**) &d_rexp, I2PI * sizeof(float));
    float* d_iexp;
    cudaMalloc((void**) &d_iexp, I2PI * sizeof(float));
    
    cudaMemcpy(d_rexp, cex_r, I2PI * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_iexp, cex_i, I2PI * sizeof(float), cudaMemcpyHostToDevice);
    
    initarraygpu<<<threadsPerGrid, threadsPerBlock>>>(d_rtcsf, d_itcsf, nnum);
    
    for(int l=0; l< __diffuse_mod_MOD_nxat; l++)
      {
	computestrucf<<<threadsPerGrid, threadsPerBlock>>>(xat[l], xat[l+__crystal_mod_MOD_cr_natoms+1], xat[l+__crystal_mod_MOD_cr_natoms+2], d_rexp, d_iexp, d_rtcsf, d_itcsf, __diffuse_mod_MOD_num[0], __diffuse_mod_MOD_num[1], __diffuse_mod_MOD_xm[0], __diffuse_mod_MOD_xm[1], __diffuse_mod_MOD_xm[2], __diffuse_mod_MOD_uin[0], __diffuse_mod_MOD_uin[1], __diffuse_mod_MOD_uin[2], __diffuse_mod_MOD_vin[0], __diffuse_mod_MOD_vin[1], __diffuse_mod_MOD_vin[2], __diffuse_mod_MOD_win[0], __diffuse_mod_MOD_win[1], __diffuse_mod_MOD_win[2]);
      }
    
    cudaMemcpy(csf_r, d_rtcsf, nnum*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(csf_i, d_itcsf, nnum*sizeof(float), cudaMemcpyDeviceToHost);
    
    
    cudaFree(d_rtcsf);
    cudaFree(d_itcsf);
    cudaFree(d_rexp);
    cudaFree(d_iexp);
    
  }
}

__global__ void computestrucf(float xat1, float xat2, float xat3, 
			      float* exp_r, float* exp_i,
			      float* tcsf_r, float* tcsf_i,
			      int num1, int num2,
			      float xm1, float xm2, float xm3,
			      float uin1, float uin2, float uin3,
			      float vin1, float vin2, float vin3, 
			      float win1, float win2, float win3)
{
  float xarg0, xincu, xincv;//, xincw;
  int iarg, iarg0, iincu, iincv;//, iincw;
  int i, j, iadd, id;
  
  xarg0 = xm1 * xat1 + xm2 * xat2 + xm3 * xat3;
  xincu = uin1 * xat1 + uin2 * xat2 + uin3 * xat3;
  xincv = vin1 * xat1 + vin2 * xat2 + vin3 * xat3;
  //xincw = win1 * xat1 + win2 * xat2 + win3 * xat3;
  iarg0 = (int)rintf(64 * I2PI * (xarg0 - (int)xarg0 + 1.));
  iincu = (int)rintf(64 * I2PI * (xincu - (int)xincu + 1.));
  iincv = (int)rintf(64 * I2PI * (xincv - (int)xincv + 1.));
  //iincw = (int)rintf(64 * I2PI * (xincw - (int)xincw + 1.));
  
  id = threadIdx.x + blockDim.x * blockIdx.x;
  if(id<num1*num2)
    {
      i = id / num1;
      j = id % num2;
      iarg = iarg0 + j * iincu + i * iincv;
      iadd = iarg >> 6;
      iadd = iadd & MASK;
      //tcsf_r[l*num1+k] += exp_r[iadd];
      //tcsf_i[l*num1+k] += exp_i[iadd];
      tcsf_r[id] += exp_r[iadd];
      tcsf_i[id] += exp_i[iadd];
    };
  __syncthreads();
}


__global__ void initarraygpu(float* array1, float* array2, int nelements)
{
  int id = threadIdx.x + blockDim.x * blockIdx.x;
  if(id<nelements)
    {
      array1[id] = 0.0;
      array2[id] = 0.0;
    };
  __syncthreads();
}

