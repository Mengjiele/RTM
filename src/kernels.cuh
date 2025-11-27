#ifndef kernels_cuh
#define kernels_cuh

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include <cuda_runtime.h>

#include "par.cuh"

// #define true    1
// #define false   0
#define EPS     0.0000000001
#define PI      3.141592653

#define BlockSize1 16	/* 1st dim block size */
#define BlockSize2 16	/* 2nd dim block size */
#define BlockSize  512	/* vector computation blocklength */
#define mm 6
#define npd     50   // boundry condition wield

#define min(x,y) x<y?x:y
#define max(x,y) x>y?x:y

#define CHECK(call) { \
    const cudaError_t err = call; \
    if (err != cudaSuccess) { \
        printf("Error: %s:%d, ", __FILE__, __LINE__); \
        printf("code:%d, reason: %s\n", err, cudaGetErrorString(err)); \
        exit(1); \
    } \
}

#define check_gpu_error(msg) { \
    cudaError_t err = cudaGetLastError(); \
    if (err != cudaSuccess) { \
        printf("Error: %s:%d, ", __FILE__, __LINE__); \
        printf("message: \"%s\", code:%d, reason: %s\n", msg, err, cudaGetErrorString(err)); \
        exit(1); \
    } \
}

__global__ void cuda_ricker_wavelet(float *wlt, float amp, float fm, float dt, int nt);
__global__ void pml_get_damp3d(float *damp, float *damp1dx, float dx, float dy, float dz, int nx, int ny, int nz, float vmin);
__global__ void cuda_add_source_lerp(bool add, float *p, float *source, float *sz, float *sx, float *sy, int nnz, int nnx, int nny);
__global__ void step_forward3d(float *p0, float *p1, float *p2, float *vv, float *damp, float dz, float dx, float dy, float dt, int nnz, int nnx, int nny);
__global__ void cuda_record_lerp(float *P, float *seis, float *gz, float *gx, float *gy, int ng, int it, int nt, int nnz, int nnx, int nny, bool record);
__global__ void cuda_cal_illum(float *illum, float *s_P, int nz, int nx, int ny);
__global__ void mute_directwave(float *seis, float *sz, const float *sx, float *sy, float *gz, float *gx, float *gy, 
								int ng, int nt, float dt, float fm, float dz, float dx, float dy, float vmute);
__global__ void save_bndr(int nnx, int nny, int nnz, int nx, int ny, int nz, int nt, float *P, float *P_bndr, bool flag);
__global__ void cross_correlate(float *mig, float *s_P, float *g_P, int nx, int ny, int nz);
__global__ void cuda_imaging(float *mig, float *illum, int nx, int ny, int nz);

__global__ void cross_correlate_ck(float *mig, float *s_P, float *g_P, int nx, int ny, int nz);
	

#ifdef __cplusplus
extern "C" {
#endif


void printPctBar(int pct, int npct, float elapsed);
void read_data(char *filename, float *tmp, int size);
void write_data(char *filename, float *tmp, int size);
void velocity_transform(float *v0, float*vv, float dt, int n1, int n2, int n3);
void read_sources(const char *filename, float *sz, float *sx, float *sy, int ns, float dz, float dx, float dy);
void read_receivers(const char *filename, float *gz, float *gx, float *gy, int ns, int is, int ng, float dz, float dx, float dy);
void window3d(float *out, float *in, int nz, int nx, int ny);
void save_wavefield(char *filename, float *temp, float *tmp, float *P, int it, int nt, int nk, int nz, int nx, int ny, int nnz, int nnx, int nny);
void load_wavefield(char *filename, float *tmp, float *g_P, float *s_P, float *mig, 
					int it, int nt, int nk, int nz, int nx, int ny, int nnz, int nnx, int nny);


void forward(float *d_vv, float *damp, float *d_wlt, float *d_sz, float *d_sx, float *d_sy, char *Recefile, 
			float dz, float dx, float dy, float dt, int nt, int ns, int ng, int nnz, int nnx, int nny);

void rtm(float *d_vv, float *damp, float *d_wlt, float *d_sz, float *d_sx, float *d_sy, char *Recefile, 
		float dz, float dx, float dy, float dt, float fm, float vmute, int nt, int ns, int ng, 
		int nz, int nx, int ny, int nnz, int nnx, int nny);
void rtm_ck(float *d_vv, float *damp, float *d_wlt, float *d_sz, float *d_sx, float *d_sy, char *Recefile, 
		float dz, float dx, float dy, float dt, float fm, float vmute, int nt, int ns, int ng, 
		int nz, int nx, int ny, int nnz, int nnx, int nny);



#ifdef __cplusplus
}
#endif









#endif