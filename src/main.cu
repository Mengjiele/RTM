//a#########################################################
//a##         	   3D Acoustic wave equation 
//a##          		forward modeling + RTM
//a##         second order constant-density equation
//a##        
//a##  (dp/dt)^2 = vp^2*[(dp/dx)^2+(dp/dy)^2+(dp/dz)^2] + f ;
//a##
//a##		boundary condition：PML
//a##		add the input acquire file
//a##                         modify by Liu Yinghui
//a##                               2025.04.18
//a#########################################################
#include "kernels.cuh"
#include "par.cuh"

//a#############################################################################################
/************************* main function for GPU Forward method ********************************/
int main( int argc, char *argv[] )
{
/*********************************variables on host*********************************************/
	int mode, nz, nx, ny, nnz, nnx, nny, nt, ns, ng;
	float *sz_h, *sx_h, *sy_h;
	float dz, dx, dy, dt, fm, amp, vmin;
	float *vp0, *vp;
/*******************************variables on devices(GPU)***************************************/
	float *d_sz, *d_sx, *d_sy;
	float *d_wlt;
	float *d_vv;
	float *damp, *damp1dx;
/***********************************************************************************************/
	cudaEvent_t start, end;
	cudaEventCreate(&start);	cudaEventCreate(&end);
	float milliseconds = 0;
	cudaEventRecord(start);
/***********************************************************************************************/
	char vpfile[250], Sourfile[250], Recefile[250];
	FILE *fp;
/***********************************************************************************************/
	parse_command_line(argc, argv);
/***********************************************************************************************/
	if(!getparint("mode", &mode))	err("must give mode= ");
	printf("================================================================================\n");
	if(mode==0) printf("####### Forward modeling \n");
	else if(mode==2) printf("#######	Reverse-time Migration \n");
/***********************************************************************************************/
	if(!getparint("nz",&nz)) err("must give nz= ");
	if(!getparint("nx",&nx)) err("must give nx= ");
	if(!getparint("ny",&ny)) err("must give ny= ");
	if(!getparfloat("dz",&dz)) err("must give dz= ");
	if(!getparfloat("dx",&dx)) err("must give dx= ");
	if(!getparfloat("dy",&dy)) err("must give dy= ");
/***********************************************************************************************/	
	if(!getparint("nt", &nt)) err("must give nt= ");
	if(!getparfloat("dt",&dt)) err("must give dt= ");
	if(!getparfloat("fm",&fm)) err("must give fm= ");
	if(!getparfloat("amp",&amp)) err("must give amp= ");
/***********************************************************************************************/	
	if(!getparint("ns", &ns)) err("must give ns= ");
	if(!getparint("ng", &ng)) err("must give ng= ");
/***********************************************************************************************/	
	if(!getparfloat("vmin",&vmin)) err("must give vmin= ");
/***********************************************************************************************/
	nnz = nz + 2*npd;
	nnx = nx + 2*npd;
	nny = ny + 2*npd;
/***********************************************************************************************/
	if (!getparstring("vpfile", vpfile)) err("must give vp_model");
	if (!getparstring("Sourfile", Sourfile)) err("must give sources");
	if (!getparstring("Recefile", Recefile)) err("must give receivers");
/********************************** allocate CPU arrays ****************************************/
	vp0					= (float*)malloc(nz*nx*ny*sizeof(float));
	vp					= (float*)malloc(nnz*nnx*nny*sizeof(float));
	
	sz_h				= (float*)malloc(ns*sizeof(float));
	sx_h				= (float*)malloc(ns*sizeof(float));
	sy_h				= (float*)malloc(ns*sizeof(float));
/********************************** read velocity file *****************************************/
	fp = fopen(vpfile, "rb");
	if(fp == NULL) {printf("###   < %s > read error!\n", vpfile);	exit(0);}
	fread(vp0, sizeof(float), nz*nx*ny, fp);
	fclose(fp);
	printf("read_vpfile done!!!!!\n");
/*************************** extension model use npd points ************************************/
	velocity_transform(vp0, vp, dt, nz, nx, ny);
/***********************************************************************************************/
	cudaSetDevice(0);// initialize device, default device=0;
	check_gpu_error("Failed to initialize device!");
/******************************* allocate GPU arrays *******************************************/
	CHECK(cudaMalloc(&d_vv,				nnz*nnx*nny*sizeof(float)));
	CHECK(cudaMalloc(&d_wlt,			nt*sizeof(float)));
	CHECK(cudaMalloc(&damp,				nnz*nnx*nny*sizeof(float)));
	CHECK(cudaMalloc(&damp1dx,			npd*sizeof(float)));
	cudaMemset(damp1dx,		0,			npd*sizeof(float));
	
	CHECK(cudaMalloc(&d_sz,				ns*sizeof(float)));
	CHECK(cudaMalloc(&d_sx,				ns*sizeof(float)));
	CHECK(cudaMalloc(&d_sy,				ns*sizeof(float)));
	check_gpu_error("Failed to allocate memory for variables!");
/***********generate nt points array of ricker with time delay and amp, fm *********************/
	cuda_ricker_wavelet<<<(nt+511)/512, 512>>>(d_wlt, amp, fm, dt, nt);
	cudaMemcpy(d_vv, vp, nnz*nnx*nny*sizeof(float), cudaMemcpyHostToDevice);
/*********************  read source point positions *******************************/
	read_sources(Sourfile, sz_h, sx_h, sy_h, ns, dz, dx, dy);
	CHECK(cudaMemcpy(d_sz, sz_h, ns*sizeof(float), cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(d_sx, sx_h, ns*sizeof(float), cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(d_sy, sy_h, ns*sizeof(float), cudaMemcpyHostToDevice));
/********************************* PML boundary coefficient ************************************/
	pml_get_damp3d<<<(nnz*nnx+511)/512,512>>>(damp, damp1dx, dx, dy, dz, nx, ny, nz, vmin);
/***********************************************************************************************/
	CHECK(cudaFree(damp1dx));
	free(vp0);		free(vp);
	free(sx_h);		free(sy_h);		free(sz_h);
/***********************************************************************************************/
	printf("--------------------------------------------------------\n");
	printf("---   \n");
/***********************************************************************************************/
	if(mode==0)		forward(d_vv, damp, d_wlt, d_sz, d_sx, d_sy, Recefile, 
							dz, dx, dy, dt, nt, ns, ng, nnz, nnx, nny);
	
	else if(mode==2)	rtm_ck(d_vv, damp, d_wlt, d_sz, d_sx, d_sy, Recefile, 
							dz, dx, dy, dt, fm, vmin, nt, ns, ng, 
							nz, nx, ny, nnz, nnx, nny);
	
	cudaEventRecord(end);
	cudaEventSynchronize(end);
	cudaEventElapsedTime(&milliseconds, start, end);
	printf("\n");
	printf("---   Complete!!!!!!!!! \n");
	printf("Time spend: %f (s)\n", milliseconds/1000.0f);
/***********************************************************************************************/	
	CHECK(cudaFree(d_vv));		CHECK(cudaFree(d_wlt));		CHECK(cudaFree(damp));
	CHECK(cudaFree(d_sx));		CHECK(cudaFree(d_sy));		CHECK(cudaFree(d_sz));

	exit (0);
}

