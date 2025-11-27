#include "kernels.cuh"
#include "par.cuh"


/************************************** CPU Function *******************************/
//a#############################################################################################
void printPctBar(int pct, int npct, float elapsed)
{
	int width = 20;  // Width of the progress bar
	float ratio = (float) (pct+1) / npct;
	int pos = width * ratio;
	
	printf("\r");
	fflush(stdout);
    
	printf("---  IS =%d/%d  ###   [", pct+1, npct);
	for (int i = 0; i < width; ++i)
	{
		if (i < pos) printf("=");
		else printf("-");
	}
	
	printf("] %d%%", (int)(ratio*100));
	printf("   ###   time spend: %f (s)", elapsed);
	fflush(stdout);
}

//a#############################################################################################
void read_data(char *filename, float *tmp, int size)
{
	FILE *fp;
	fp = fopen(filename, "rb");
	if(fp==NULL) { fprintf(stderr,"read_data, error opening file\n"); exit(1);}
	fread(tmp, sizeof(float), size, fp);
	fclose(fp);
}

//a#############################################################################################
void write_data(char *filename, float *tmp, int size)
{
	FILE *fp;
	fp = fopen(filename, "wb");
	if(fp==NULL) { fprintf(stderr,"write_data, error opening file\n"); exit(1);}
	fwrite(tmp, sizeof(float), size, fp);
	fclose(fp);
}

//a#############################################################################################
void velocity_transform(float *v0, float*vv, float dt, int n1, int n2, int n3)
 /*< velocit2 transform: vv=v0*dt; vv<--vv^2 >*/
{
  int i1, i2, i3, nn1, nn2, nn3;
  float tmp;

  nn1=n1+2*npd;
  nn2=n2+2*npd;
  nn3=n3+2*npd;

  // inner zone
  for(i3=0; i3<n3; i3++){//y
    for(i2=0; i2<n2; i2++){//x
      for(i1=0; i1<n1; i1++){//z
	tmp=v0[i1+n1*i2+n1*n2*i3];//*dt;
	vv[(i1+npd)+nn1*(i2+npd)+nn1*nn2*(i3+npd)]=tmp;//*tmp;
      }
    }
  }  
    //top & down 
    for(i3=0; i3<nn3; i3++){//y
	for(i2=0; i2<nn2; i2++){//x
	    for (i1=0; i1<npd; i1++){//z
		vv[i1+nn1*i2+nn1*nn2*i3]=vv[npd+nn1*i2+nn1*nn2*i3];
		vv[(nn1-i1-1)+nn1*i2+nn1*nn2*i3]=vv[(nn1-npd-1)+nn1*i2+nn1*nn2*i3];
	    }
	}
    }

    //left & right
    for(i3=0; i3<nn3; i3++){//y
	for(i2=0; i2<npd; i2++){//x
	    for (i1=0; i1<nn1; i1++){//z
		vv[i1+nn1*i2+nn1*nn2*i3]=vv[i1+nn1*(npd)+nn1*nn2*i3];
		vv[i1+nn1*(nn2-i2-1)+nn1*nn2*i3]=vv[i1+nn1*(nn2-npd-1)+nn1*nn2*i3];
	    }
	}
    }
    //front & back
    for(i3=0; i3<npd; i3++){//y
	for(i2=0; i2<nn2; i2++){//x
	    for(i1=0; i1<nn1; i1++){//z
		vv[i1+nn1*i2+nn1*nn2*i3]=vv[i1+nn1*i2+nn1*nn2*(npd)];
		vv[i1+nn1*i2+nn1*nn2*(nn3-1-i3)]=vv[i1+nn1*i2+nn1*nn2*(nn3-npd-1)];
	    }
	}
    }
}

//a#############################################################################################
void read_sources(const char *filename, float *sz, float *sx, float *sy, int ns, float dz, float dx, float dy)
{
	FILE *fp = fopen(filename, "r"); 
	if (!fp) { perror("read_sources fopen"); exit(1); }
	
	int tmp;
	float z, x, y;
	
	char head[256];
	if(fgets(head, sizeof(head), fp) == NULL)//skip a line at the beginning of the receivers file
	{
		fprintf(stderr, "Error: File is empty or cannot read head.\n");
		fclose(fp);
		exit(1);
	}
	
	for(int i=0; i<ns; i++)
	{
		fscanf(fp, "%d %f %f %f", &tmp, &z, &x, &y);
		sz[i] = z/dz;
		sx[i] = x/dx;
		sy[i] = y/dy;
	}
	fclose(fp);
}

//a#############################################################################################
void read_receivers(const char *filename, float *gz, float *gx, float *gy, int ns, int is, int ng, float dz, float dx, float dy)
{
	FILE *fp = fopen(filename, "r");
	if (!fp) { perror("read_receivers fopen"); exit(1); }
	
	int tmp;
	float z, x, y;
	
	char head[256];
	if(fgets(head, sizeof(head), fp) == NULL)//skip a line at the beginning of the receivers file
	{
		fprintf(stderr, "Error: File is empty or cannot read head.\n");
		fclose(fp);
		exit(1);
	}
	
	for (int i = 0; i < is * ng; i++)//skip the data of the (is*ng)
	{
		if (fscanf(fp, "%*d %*f %*f %*f") == EOF){
			if (feof(fp)){
				fprintf(stderr, "Error: File ended early. Expected %d lines, but got only %d.\n", is*ng, i);
			}else{
				perror("Error reading receiver data");
			}
			fclose(fp);
			exit(1);
		}
	}
	
	for(int i=0; i<ng; i++)
	{
		fscanf(fp, "%d %f %f %f", &tmp, &z, &x, &y);
		gz[i] = z/dz;
		gx[i] = x/dx;
		gy[i] = y/dy;
	}
	fclose(fp);
}

//a#############################################################################################
void window3d(float *out, float *in, int nz, int nx, int ny)
/*< window a 3d subvolume >*/
{
	int iz, ix, iy;
	int nnz=nz+2*npd;//z
	int nnx=nx+2*npd;//x
	
	for(iy=0; iy<ny; iy++)
	for(ix=0; ix<nx; ix++)
	for(iz=0; iz<nz; iz++)
	{
		out[iz+nz*ix+nz*nx*iy] = in[(iz+npd)+nnz*(ix+npd)+nnz*nnx*(iy+npd)];
	}
}

/************************************** GPU Function *******************************/
//a#############################################################################################
__global__ void cuda_ricker_wavelet(float *wlt, float amp, float fm, float dt, int nt)
/*< generate ricker wavelet with time deley >*/
{
	int it=threadIdx.x+blockDim.x*blockIdx.x;
	if (it<nt){
		float tmp = PI*fm*fabsf(it*dt-1.0/fm);//delay the wavelet to exhibit all waveform
		tmp *=tmp;
		wlt[it]= amp*(1.0-2.0*tmp)*expf(-tmp);// ricker wavelet at time: t=nt*dt
	}
}

//a#############################################################################################
__global__ void pml_get_damp3d(float *damp, float *damp1dx, float dx, float dy, float dz, int nx, int ny, int nz, float vmin)
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	
	int nnz = nz+2*npd;
	int nnx = nx+2*npd;
	int nny = ny+2*npd;
	
	int iy=id/(nnz*nnx);
	int ix=(id%(nnz*nnx))/nnz;
	int iz=(id%(nnz*nnx))%nnz;
	
	float xa;
	float a = (npd-1)*dx;
	float kappa = 3.0*vmin*log(100000.0)/(2.0*a);
	
	for(int i=0; i<npd; i++)
	{
		xa=(float)i*dx/a;
		damp1dx[i]=kappa*xa*xa;
	}
	for(iy=0; iy<nny; iy++)
	{
		if(ix<npd && iz<nnz)
		{
			damp[iz+(npd-ix-1)*nnz+iy*nnz*nnx]=damp1dx[ix];
			damp[iz+(nx+npd+ix)*nnz+iy*nnz*nnx]=damp1dx[ix];
		}
	}
	for(iy=0; iy<nny; iy++)
	{
		for(iz=0; iz<npd; iz++)
		{
			for(ix=npd-iz; ix<nx+npd+iz; ix++)
			{
				damp[(npd-iz-1)+ix*nnz+iy*nnz*nnx]=damp1dx[iz];
				damp[(nz+npd+iz)+ix*nnz+iy*nnz*nnx]=damp1dx[iz];
			}
		}
	}
	for(iy=0; iy<npd; iy++)
	{
		for(ix=npd-iy; ix<nx+npd+iy; ix++)
		{
			for(iz=npd-iy; iz<nz+npd+iy; iz++)
			{
				damp[iz+ix*nnz+(npd-iy-1)*nnz*nnx]=damp1dx[iy];
				damp[iz+ix*nnz+(ny+npd+iy)*nnz*nnx]=damp1dx[iy];
			}
		}
	}
}

//a#############################################################################################
__global__ void cuda_add_source_lerp(bool add, float *p, float *source, float *sz, float *sx, float *sy, int nnz, int nnx, int nny)
{
	float z = sz[0] + npd;
	float x = sx[0] + npd;
	float y = sy[0] + npd;
	
	int z0 = min(max((int)z, 0), nnz-2);
	int x0 = min(max((int)x, 0), nnx-2);
	int y0 = min(max((int)y, 0), nny-2);
	int z1 = z0 + 1;
	int x1 = x0 + 1;
	int y1 = y0 + 1;
	
	float wz = min(max(z - z0, 0.0f), 1.0f);
	float wx = min(max(x - x0, 0.0f), 1.0f);
	float wy = min(max(y - y0, 0.0f), 1.0f);
	
	float w000 = (1-wz)*(1-wx)*(1-wy);
	float w100 = wz*(1-wx)*(1-wy);
	float w010 = (1-wz)*wx*(1-wy);
	float w110 = wz*wx*(1-wy);
	float w001 = (1-wz)*(1-wx)*wy;
	float w101 = wz*(1-wx)*wy;
	float w011 = (1-wz)*wx*wy;
	float w111 = wz*wx*wy;
	
	int idx000 = z0 + nnz*(x0 + nnx*y0);
	int idx100 = z1 + nnz*(x0 + nnx*y0);
	int idx010 = z0 + nnz*(x1 + nnx*y0);
	int idx110 = z1 + nnz*(x1 + nnx*y0);
	int idx001 = z0 + nnz*(x0 + nnx*y1);
	int idx101 = z1 + nnz*(x0 + nnx*y1);
	int idx011 = z0 + nnz*(x1 + nnx*y1);
	int idx111 = z1 + nnz*(x1 + nnx*y1);
	
	float src_val = source[0];
	
	if (add)
	{
		atomicAdd(&p[idx000], src_val * w000);
		atomicAdd(&p[idx100], src_val * w100);
		atomicAdd(&p[idx010], src_val * w010);
		atomicAdd(&p[idx110], src_val * w110);
		atomicAdd(&p[idx001], src_val * w001);
		atomicAdd(&p[idx101], src_val * w101);
		atomicAdd(&p[idx011], src_val * w011);
		atomicAdd(&p[idx111], src_val * w111);
    }
	else
	{
        atomicAdd(&p[idx000], -src_val * w000);
		atomicAdd(&p[idx100], -src_val * w100);
		atomicAdd(&p[idx010], -src_val * w010);
		atomicAdd(&p[idx110], -src_val * w110);
		atomicAdd(&p[idx001], -src_val * w001);
		atomicAdd(&p[idx101], -src_val * w101);
		atomicAdd(&p[idx011], -src_val * w011);
		atomicAdd(&p[idx111], -src_val * w111);
	}
	
}


//a#############################################################################################
__global__ void step_forward3d(float *p0, float *p1, float *p2, float *vv, float *damp, float dz, float dx, float dy, float dt, int nnz, int nnx, int nny)
{
	int iz = blockIdx.x * blockDim.x + threadIdx.x;
	int ix = blockIdx.y * blockDim.y + threadIdx.y;
	int iy = blockIdx.z * blockDim.z + threadIdx.z;
	
	if (ix >= nnx || iy >= nny || iz >= nnz)
		return;
	
	int id = iz + ix * nnz + iy * nnz * nnx;
	
	if (ix < mm || ix >= nnx - mm ||
		iy < mm || iy >= nny - mm ||
		iz < mm || iz >= nnz - mm)
		return;
		
	float dtx = dt / dx;
	float dty = dt / dy;
	float dtz = dt / dz;
	
	float alpha_x = vv[id]*vv[id] * dtx * dtx;
	float alpha_y = vv[id]*vv[id] * dty * dty;
	float alpha_z = vv[id]*vv[id] * dtz * dtz;
	
	float kappa = damp[id] * dt;
	float a1 = 2.0f - kappa * kappa;
	float a2 = 1.0f - kappa;
	float a3 = 1.0f + kappa;
	
	const float c0 = -2.98278f,		c1 = 1.71429f,		c2 = -0.267857f,		c3 = 0.0529101f;
	const float c4 = -0.00892857f,	c5 = 0.00103896f,	c6 = -0.0000601251f;
	
	float dz_term =	c0 * p1[id]	+
					c1 * (p1[id+1] + p1[id-1]) +
					c2 * (p1[id+2] + p1[id-2]) +
					c3 * (p1[id+3] + p1[id-3]) +
					c4 * (p1[id+4] + p1[id-4]) +
					c5 * (p1[id+5] + p1[id-5]) +
					c6 * (p1[id+6] + p1[id-6]);
	
	float dx_term =	c0 * p1[id] +
					c1 * (p1[id+nnz] + p1[id-nnz]) +
					c2 * (p1[id+2*nnz] + p1[id-2*nnz]) +
					c3 * (p1[id+3*nnz] + p1[id-3*nnz]) +
					c4 * (p1[id+4*nnz] + p1[id-4*nnz]) +
					c5 * (p1[id+5*nnz] + p1[id-5*nnz]) +
					c6 * (p1[id+6*nnz] + p1[id-6*nnz]);
	
	int nnzx = nnz * nnx;
	
	float dy_term =	c0 * p1[id] +
					c1 * (p1[id+nnzx] + p1[id-nnzx]) +
					c2 * (p1[id+2*nnzx] + p1[id-2*nnzx]) +
					c3 * (p1[id+3*nnzx] + p1[id-3*nnzx]) +
					c4 * (p1[id+4*nnzx] + p1[id-4*nnzx]) +
					c5 * (p1[id+5*nnzx] + p1[id-5*nnzx]) +
					c6 * (p1[id+6*nnzx] + p1[id-6*nnzx]);
					
	p2[id] = (a1*p1[id] - a2*p0[id] + alpha_z*dz_term + alpha_x*dx_term + alpha_y*dy_term)/a3;
	
}

//a#############################################################################################
__global__ void cuda_record_lerp(float *P, float *seis, float *gz, float *gx, float *gy, int ng, int it, int nt, int nnz, int nnx, int nny, bool record)
{
	int id = threadIdx.x + blockDim.x * blockIdx.x;
	
	
	if (id<ng)
	{
		float z = gz[id] + npd;
		float x = gx[id] + npd;
		float y = gy[id] + npd;
		
		int z0 = min(max((int)z, 0), nnz-2);
		int x0 = min(max((int)x, 0), nnx-2);
		int y0 = min(max((int)y, 0), nny-2);
		int z1 = z0 + 1;
		int x1 = x0 + 1;
		int y1 = y0 + 1;
		
		float wz = min(max(z - z0, 0.0f), 1.0f);
		float wx = min(max(x - x0, 0.0f), 1.0f);
		float wy = min(max(y - y0, 0.0f), 1.0f);
		
		float v000 = P[z0 + nnz*(x0 + nnx*y0)];
		float v100 = P[z1 + nnz*(x0 + nnx*y0)];
		float v010 = P[z0 + nnz*(x1 + nnx*y0)];
		float v110 = P[z1 + nnz*(x1 + nnx*y0)];
		float v001 = P[z0 + nnz*(x0 + nnx*y1)];
		float v101 = P[z1 + nnz*(x0 + nnx*y1)];
		float v011 = P[z0 + nnz*(x1 + nnx*y1)];
		float v111 = P[z1 + nnz*(x1 + nnx*y1)];
		
		float interpolated = 
			v000 * (1-wz)*(1-wx)*(1-wy) +
			v100 * wz*(1-wx)*(1-wy) +
			v010 * (1-wz)*wx*(1-wy) +
			v110 * wz*wx*(1-wy) +
			v001 * (1-wz)*(1-wx)*wy +
			v101 * wz*(1-wx)*wy +
			v011 * (1-wz)*wx*wy +
			v111 * wz*wx*wy;
			
		if (record)
		{
			seis[it+id*nt] = interpolated;
		}
		else
		{
			float val = seis[it + id*nt];
			atomicAdd(&P[z0 + nnz*(x0 + nnx*y0)], val * (1-wz)*(1-wx)*(1-wy));
			atomicAdd(&P[z1 + nnz*(x0 + nnx*y0)], val * wz*(1-wx)*(1-wy));
			atomicAdd(&P[z0 + nnz*(x1 + nnx*y0)], val * (1-wz)*wx*(1-wy));
			atomicAdd(&P[z1 + nnz*(x1 + nnx*y0)], val * wz*wx*(1-wy));
			atomicAdd(&P[z0 + nnz*(x0 + nnx*y1)], val * (1-wz)*(1-wx)*wy);
			atomicAdd(&P[z1 + nnz*(x0 + nnx*y1)], val * wz*(1-wx)*wy);
			atomicAdd(&P[z0 + nnz*(x1 + nnx*y1)], val * (1-wz)*wx*wy);
			atomicAdd(&P[z1 + nnz*(x1 + nnx*y1)], val * wz*wx*wy);
		}
		
	}
}

//a#############################################################################################
__global__ void cuda_cal_illum(float *illum, float *s_P, int nz, int nx, int ny)
/*< calculate the source lighting matrix >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	
	int iy = id/(nz*nx);
	int ix = (id%(nz*nx))/nz;
	int iz = (id%(nz*nx))%nz;
	
	int nnz = nz + 2*npd;
	int nnx = nx + 2*npd;
	
	if(id<nz*nx*ny)
	{
		if(ix<nx && iy<ny && iz<nz)
		{
			illum[id]+=s_P[(iz+npd)+(ix+npd)*nnz+(iy+npd)*nnz*nnx]*s_P[(iz+npd)+(ix+npd)*nnz+(iy+npd)*nnz*nnx];
		}
	}
}

__global__ void mute_directwave(float *seis, float *sz, const float *sx, float *sy, float *gz, float *gx, float *gy, 
								int ng, int nt, float dt, float fm, float dz, float dx, float dy, float vmute)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if (id >= ng) return;
	
	float mu_z = (*sz) - gz[id];
	float mu_x = (*sx) - gx[id];
	float mu_y = (*sy) - gy[id];
	float dist = sqrtf(mu_z*mu_z*dz*dz + mu_x*mu_x*dx*dx + mu_y*mu_y*dy*dy);
	
	float t_mute = dist/vmute;
	int it_mute = (int)(t_mute/dt) + (int)(2.0f/(dt*fm));//(int)(1.8f/(dt*fm))
//	int it_taper = (int)(0.05f / dt);//it_taper = (int)(taper_time / dt);
	
	for (int it = 0; it < nt; ++it)
	{
		int idx = it + id * nt;
		if (it < it_mute)
			seis[idx] = 0.0f;
	/*	else if (it < it_mute + it_taper)
		{
		//	float taper = (float)(it - it_mute) / it_taper;
			float taper = 0.5f * (1.0f + cosf(PI*(float)(it-it_mute)/(float)it_taper));
			seis[idx] *= taper;
		}*/
	}
	
}

//a#############################################################################################
__global__ void save_bndr(int nnx, int nny, int nnz, int nx, int ny, int nz, int nt, float *P, float *P_bndr, bool flag)
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int ix,iy,iz;
	
	if(id<2*nx*ny+2*nz*ny+2*nx*nz)
	{
		if(flag)/////////////////////////////////save boundary
		{
			if(id<nx*ny){//up
				ix=id%nx;
				iy=id/nx;
				P_bndr[id]=P[npd-1+nnz*(ix+npd)+nnz*nnx*(iy+npd)];
				
			}else if(id>=nx*ny&&id<(2*nx*ny)){//down
				ix=(id-nx*ny)%nx;
				iy=(id-nx*ny)/nx;
				P_bndr[id]=P[npd+nz+nnz*(ix+npd)+nnz*nnx*(iy+npd)];
			
			}else if(id>=(2*nx*ny)&&id<(2*nx*ny+nz*ny)){//left
				iz=(id-2*nx*ny)%nz;
				iy=(id-2*nx*ny)/nz;
				P_bndr[id]=P[npd+iz+nnz*(npd-1)+nnz*nnx*(iy+npd)];

			}else if(id>=(2*nx*ny+nz*ny)&&id<(2*nx*ny+2*nz*ny)){//right
				iz=(id-2*nx*ny-nz*ny)%nz;
				iy=(id-2*nx*ny-nz*ny)/nz;
				P_bndr[id]=P[npd+iz+nnz*(nx+npd)+nnz*nnx*(iy+npd)];

			}else if(id>=(2*nx*ny+2*nz*ny)&&id<(2*nx*ny+2*nz*ny+nx*nz)){//front
				iz=(id-2*nx*ny-2*nz*ny)%nz;
				ix=(id-2*nx*ny-2*nz*ny)/nz;
				P_bndr[id]=P[npd+iz+nnz*(ix+npd)+nnz*nnx*(npd-1)];

			}else if(id>=(2*nx*ny+2*nz*ny+nx*nz)&&id<(2*nx*ny+2*nz*ny+2*nx*nz)){//back
				iz=(id-2*nx*ny-2*nz*ny-nx*nz)%nz;
				ix=(id-2*nx*ny-2*nz*ny-nx*nz)/nz;
				P_bndr[id]=P[npd+iz+nnz*(ix+npd)+nnz*nnx*(npd+ny)];
			}
		
		}else{
			if(id<nx*ny){//up
				ix=id%nx;
				iy=id/nx;
				P[npd-1+nnz*(ix+npd)+nnz*nnx*(iy+npd)]=P_bndr[id];
			
			}else if(id>=nx*ny&&id<(2*nx*ny)){//down
				ix=(id-nx*ny)%nx;
				iy=(id-nx*ny)/nx;
				P[npd+nz+nnz*(ix+npd)+nnz*nnx*(iy+npd)]=P_bndr[id];

			}else if(id>=(2*nx*ny)&&id<(2*nx*ny+nz*ny)){//left
				iz=(id-2*nx*ny)%nz;
				iy=(id-2*nx*ny)/nz;
				P[npd+iz+nnz*(npd-1)+nnz*nnx*(iy+npd)]=P_bndr[id];

			}else if(id>=(2*nx*ny+nz*ny)&&id<(2*nx*ny+2*nz*ny)){//right
				iz=(id-2*nx*ny-nz*ny)%nz;
				iy=(id-2*nx*ny-nz*ny)/nz;
				P[npd+iz+nnz*(nx+npd)+nnz*nnx*(iy+npd)]=P_bndr[id];

			}else if(id>=(2*nx*ny+2*nz*ny)&&id<(2*nx*ny+2*nz*ny+nx*nz)){//front
				iz=(id-2*nx*ny-2*nz*ny)%nz;
				ix=(id-2*nx*ny-2*nz*ny)/nz;
				P[npd+iz+nnz*(ix+npd)+nnz*nnx*(npd-1)]=P_bndr[id];

			}else if(id>=(2*nx*ny+2*nz*ny+nx*nz)&&id<(2*nx*ny+2*nz*ny+2*nx*nz)){//back
				iz=(id-2*nx*ny-2*nz*ny-nx*nz)%nz;
				ix=(id-2*nx*ny-2*nz*ny-nx*nz)/nz;
				P[npd+iz+nnz*(ix+npd)+nnz*nnx*(npd+ny)]=P_bndr[id];
			}
		}
	}       
}

//a#############################################################################################
__global__ void cross_correlate(float *mig, float *s_P, float *g_P, int nx, int ny, int nz)
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	
	int iy = id/(nz*nx);
	int ix = (id%(nz*nx))/nz;
	int iz = (id%(nz*nx))%nz;
	
	int nnz = nz + 2*npd;
	int nnx = nx + 2*npd;
	
	if(id<nz*nx*ny)
	{
		if(ix<nx && iy<ny &&iz<nz)
		{
			mig[id]+=s_P[(iz+npd)+(ix+npd)*nnz+(iy+npd)*nnz*nnx]*g_P[(iz+npd)+(ix+npd)*nnz+(iy+npd)*nnz*nnx];
		}
	}
}

//a#############################################################################################
__global__ void cross_correlate_ck(float *mig, float *s_P, float *g_P, int nx, int ny, int nz)
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	
	int iy = id/(nz*nx);
	int ix = (id%(nz*nx))/nz;
	int iz = (id%(nz*nx))%nz;
	
	int nnz = nz + 2*npd;
	int nnx = nx + 2*npd;
	
	if(id<nz*nx*ny)
	{
		if(ix<nx && iy<ny &&iz<nz)
		{
			mig[id]+=s_P[iz+ix*nz+iy*nz*nx]*g_P[(iz+npd)+(ix+npd)*nnz+(iy+npd)*nnz*nnx];
		}
	}
}

//a#############################################################################################
__global__ void cuda_imaging(float *mig, float *illum, int nx, int ny, int nz)
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	
//	int iy = id/(nz*nx);
//	int ix = (id%(nz*nx))/nz;
//	int iz = (id%(nz*nx))%nz;
	
	if(id<nz*nx*ny)
	{
		if(illum[id]!=0)
		{
			mig[id]/=illum[id];
		}
	}
	
}

//a#############################################################################################
void save_wavefield(char *filename, float *temp, float *tmp, float *P, int it, int nt, int nk, 
					int nz, int nx, int ny, int nnz, int nnx, int nny)
{
	if((it%nk==0) ||(it==nt-1))
	{
		CHECK(cudaMemcpy(temp, P, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost));
		window3d(tmp, temp, nz, nx, ny);
		sprintf(filename, "./data/Output/Wavefield/check_wavefield_%d.dat", it);
		write_data(filename, tmp, nz*nx*ny);
	}
}

//a#############################################################################################
void load_wavefield(char *filename, float *tmp, float *g_P, float *s_P, float *mig, 
					int it, int nt, int nk, int nz, int nx, int ny, int nnz, int nnx, int nny)
{
	if((it%nk==0) ||(it==nt-1))
	{	
		sprintf(filename, "./data/Output/Wavefield/check_wavefield_%d.dat", it);
		read_data(filename, tmp, nz*nx*ny);
		CHECK(cudaMemcpy(s_P, tmp, nz*nx*ny*sizeof(float), cudaMemcpyHostToDevice));
		cross_correlate_ck<<<(nz*nx*ny+511)/512,512>>>(mig, s_P, g_P, nx, ny, nz);
	}
}

//a#############################################################################################
void forward(float *d_vv, float *damp, float *d_wlt, float *d_sz, float *d_sx, float *d_sy, char *Recefile, 
			float dz, float dx, float dy, float dt, int nt, int ns, int ng, int nnz, int nnx, int nny)
{
	int it, is;
	float *gz_h, *gx_h, *gy_h;
	float *p_IO, *temp;
	float *p_cal;
	float *d_gz, *d_gx, *d_gy;
	float *d_p0, *d_p1, *d_p2, *ptr;
	
	cudaEvent_t is_t0, is_t1, ns_t0, ns_t1;
	cudaEventCreate(&is_t0);	cudaEventCreate(&is_t1);
	cudaEventCreate(&ns_t0);	cudaEventCreate(&ns_t1);
	float seconds = 0;
	float is_elapsed;
	
	char filename[250];
	FILE *fp;
	
	temp				= (float*)malloc(nnz*nnx*nny*sizeof(float));
	p_IO				= (float*)malloc(nt*ng*sizeof(float));
	
	gz_h				= (float*)malloc(ng*sizeof(float));
	gx_h				= (float*)malloc(ng*sizeof(float));
	gy_h				= (float*)malloc(ng*sizeof(float));
/******************************* allocate GPU arrays *******************************************/
	cudaSetDevice(0);// initialize device, default device=0;
	check_gpu_error("Failed to initialize device!");
	
	dim3 dimb(8, 8, 8);
	dim3 dimg((nnz+7)/8, (nnx+7)/8, (nny+7)/8);
//*****< forward & backward & receivers wavefield >
	CHECK(cudaMalloc(&d_p0,			nnz*nnx*nny*sizeof(float)));
	CHECK(cudaMalloc(&d_p1,			nnz*nnx*nny*sizeof(float)));
	CHECK(cudaMalloc(&d_p2,			nnz*nnx*nny*sizeof(float)));
	
	CHECK(cudaMalloc(&p_cal,		nt*ng*sizeof(float)));
	
	CHECK(cudaMalloc(&d_gz,			ng*sizeof(float)));
	CHECK(cudaMalloc(&d_gx,			ng*sizeof(float)));
	CHECK(cudaMalloc(&d_gy,			ng*sizeof(float)));
	check_gpu_error("Failed to allocate memory for variables!");
/***********************************************************************************************/
	cudaEventRecord(ns_t0);
	printf("===   Forward modeling is start\n");
/**************************** begin ns loop  *************************************/
	for(is=0; is<300; is++)
	{
		cudaEventRecord(is_t0);
	/****************** initialize arrays with 0 *************************************/
		CHECK(cudaMemset(d_p0,	0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(d_p1,	0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(d_p2,	0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(p_cal,	0,		nt*ng*sizeof(float)));
	/*********************  set  receivers point positions ***************************/
		read_receivers(Recefile, gz_h, gx_h, gy_h, ns, is, ng, dz, dx, dy);
		CHECK(cudaMemcpy(d_gz, gz_h, ng*sizeof(float), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_gx, gx_h, ng*sizeof(float), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_gy, gy_h, ng*sizeof(float), cudaMemcpyHostToDevice));
	/****************************** Forward *************************************/
		for(it=0; it<nt; it++)
		{
			cuda_add_source_lerp<<<1,1>>>(true, d_p1, &d_wlt[it], &d_sz[is], &d_sx[is], &d_sy[is], nnz, nnx, nny);
			step_forward3d<<<dimg,dimb>>>(d_p0, d_p1, d_p2, d_vv, damp, dz, dx, dy, dt, nnz, nnx, nny);
			ptr=d_p0;	d_p0=d_p1; d_p1=d_p2; d_p2=ptr;
			cuda_record_lerp<<<(ng+511)/512, 512>>>(d_p1, p_cal, d_gz, d_gx, d_gy, ng, it, nt, nnz, nnx, nny, true);
			
			if(is==0&&it==500)
			{
				sprintf(filename,"./data/Output/snap.dat");
				fp=fopen(filename,"wb");
				cudaMemcpy(temp, d_p1, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
				fwrite(&temp[0], sizeof(float), nnz*nnx*nny, fp);
				fclose(fp);
			}
			
		}//it loop end
		
		cudaMemcpy(p_IO, p_cal, nt*ng*sizeof(float), cudaMemcpyDeviceToHost);
		sprintf(filename, "./data/Input/Shot/%dshot_obs.dat", is+1);
		write_data(filename, p_IO, nt*ng);
		
		cudaEventRecord(is_t1);
		cudaEventSynchronize(is_t1);
		cudaEventElapsedTime(&seconds, is_t0, is_t1);
		is_elapsed = seconds/1000.0f;
		printPctBar(is, ns, is_elapsed);
		usleep(50000);// Sleep for 50 milliseconds
	}//is loop end
	cudaEventRecord(ns_t1);
	cudaEventSynchronize(ns_t1);
	cudaEventElapsedTime(&seconds, ns_t0, ns_t1);
	printf("\n");
	printf("total %d shots: %f (s)\n", ns, seconds/1000.0f);
	printf("---   The forward is over    \n");
	
	free(gx_h);		free(gy_h);		free(gz_h);
	free(p_IO);		free(temp);
	CHECK(cudaFree(d_p0));		CHECK(cudaFree(d_p1));		CHECK(cudaFree(d_p2));
	CHECK(cudaFree(p_cal));
	CHECK(cudaFree(d_gx));		CHECK(cudaFree(d_gy));		CHECK(cudaFree(d_gz));
}

//a#############################################################################################
void rtm(float *d_vv, float *damp, float *d_wlt, float *d_sz, float *d_sx, float *d_sy, char *Recefile, 
		float dz, float dx, float dy, float dt, float fm, float vmute, int nt, int ns, int ng, 
		int nz, int nx, int ny, int nnz, int nnx, int nny)
{
	int it, is;
	float *gz_h, *gx_h, *gy_h;
	float *p_IO, *tmp, *temp, *ptr;
	float *p_cal;
	float *d_gz, *d_gx, *d_gy;
	float *s_P0, *s_P1, *s_P2;
	float *g_P0, *g_P1, *g_P2;
	float *h_bndr, *d_bndr;
	float *illum, *mig;
	
	
	cudaEvent_t is_t0, is_t1, ns_t0, ns_t1;
	cudaEventCreate(&is_t0);	cudaEventCreate(&is_t1);
	cudaEventCreate(&ns_t0);	cudaEventCreate(&ns_t1);
	float seconds = 0;
	float is_elapsed;
	
	char filename[250];
	FILE *fp;
/***********************************************************************************************/
	int nb = 2*nz*nx + 2*nz*ny + 2*nx*ny;
/********************************** allocate CPU arrays ****************************************/
	tmp					= (float*)malloc(nz*nx*ny*sizeof(float));
	temp				= (float*)malloc(nnz*nnx*nny*sizeof(float));
	p_IO				= (float*)malloc(nt*ng*sizeof(float));
	h_bndr				= (float*)malloc(nb*sizeof(float));
	
	gz_h				= (float*)malloc(ng*sizeof(float));
	gx_h				= (float*)malloc(ng*sizeof(float));
	gy_h				= (float*)malloc(ng*sizeof(float));
/******************************* allocate GPU arrays *******************************************/
	cudaSetDevice(0);// initialize device, default device=0;
	check_gpu_error("Failed to initialize device!");
	
	dim3 dimb(8, 8, 8);
	dim3 dimg((nnz+7)/8, (nnx+7)/8, (nny+7)/8);
/******************************* allocate GPU arrays *******************************************/
//*****< forward & backward & receivers wavefield >
	CHECK(cudaMalloc(&s_P0,			nnz*nnx*nny*sizeof(float)));
	CHECK(cudaMalloc(&s_P1,			nnz*nnx*nny*sizeof(float)));
	CHECK(cudaMalloc(&s_P2,			nnz*nnx*nny*sizeof(float)));
	CHECK(cudaMalloc(&g_P0,			nnz*nnx*nny*sizeof(float)));
	CHECK(cudaMalloc(&g_P1,			nnz*nnx*nny*sizeof(float)));
	CHECK(cudaMalloc(&g_P2,			nnz*nnx*nny*sizeof(float)));
	
	CHECK(cudaMalloc(&p_cal,		nt*ng*sizeof(float)));
	
	CHECK(cudaMalloc(&d_gz,			ng*sizeof(float)));
	CHECK(cudaMalloc(&d_gx,			ng*sizeof(float)));
	CHECK(cudaMalloc(&d_gy,			ng*sizeof(float)));
	
	CHECK(cudaMalloc(&d_bndr,		nb*sizeof(float)));
//*****< The is & ns gradient ,lighting matrix >
	CHECK(cudaMalloc(&illum,		nz*nx*ny*sizeof(float)));
	CHECK(cudaMalloc(&mig,			nz*nx*ny*sizeof(float)));
/***********************************************************************************************/
	CHECK(cudaMemset(illum,		0,		nz*nx*ny*sizeof(float)));
	CHECK(cudaMemset(mig,		0,		nz*nx*ny*sizeof(float)));
	printf("===   Reverse-time Migration is start\n");
	cudaEventRecord(ns_t0);
/**************************** begin ns loop  *************************************/
	for(is=0; is<ns; is++)
	{
		cudaEventRecord(is_t0);
		CHECK(cudaMemset(s_P0,		0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(s_P1,		0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(s_P2,		0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(p_cal,		0,		nt*ng*sizeof(float)));
		CHECK(cudaMemset(d_bndr,	0,		nb*sizeof(float)));
	/*********************  set  receivers point positions ***************************/
		read_receivers(Recefile, gz_h, gx_h, gy_h, ns, is, ng, dz, dx, dy);
		CHECK(cudaMemcpy(d_gz, gz_h, ng*sizeof(float), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_gx, gx_h, ng*sizeof(float), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_gy, gy_h, ng*sizeof(float), cudaMemcpyHostToDevice));
	/****************************** Forward *************************************/
		FILE *fp_bndr;
		sprintf(filename, "./data/Output/bndr.dat");
		fp_bndr = fopen(filename, "wb");
		for(it=0; it<nt; it++)
		{
			cuda_add_source_lerp<<<1,1>>>(true, s_P1, &d_wlt[it], &d_sz[is], &d_sx[is], &d_sy[is], nnz, nnx, nny);
			step_forward3d<<<dimg,dimb>>>(s_P0, s_P1, s_P2, d_vv, damp, dz, dx, dy, dt, nnz, nnx, nny);
			ptr=s_P0;	s_P0=s_P1;	s_P1=s_P2;	s_P2=ptr;
		/**************** write boundary wavefield to arrays ********************/
			save_bndr<<<(nb+511)/512, 512>>>(nnx, nny, nnz, nx, ny, nz, nt, s_P0, d_bndr, true);
			CHECK(cudaMemcpy(h_bndr, d_bndr, nb*sizeof(float), cudaMemcpyDeviceToHost));
			fwrite(h_bndr, sizeof(float), nb, fp_bndr);
		/************************** illumination ********************************/
			cuda_cal_illum<<<(nz*nx*ny+511)/512,512>>>(illum, s_P1, nz, nx, ny);
			cuda_record_lerp<<<(ng+511)/512, 512>>>(s_P1, p_cal, d_gz, d_gx, d_gy, ng, it, nt, nnz, nnx, nny, true);
			
			if(is==0&&it==1000)
			{
				sprintf(filename,"./data/Output/snapf.dat");
				fp=fopen(filename,"wb");
				cudaMemcpy(temp, s_P1, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
				fwrite(&temp[0], sizeof(float), nnz*nnx*nny, fp);
				fclose(fp);
			}		
		}//it loop end
		fclose(fp_bndr);
		
		sprintf(filename, "./data/Input/Shot/%dshot_obs.dat", is+1);
		read_data(filename, p_IO, nt*ng);
		CHECK(cudaMemcpy(p_cal, p_IO, nt*ng*sizeof(float), cudaMemcpyHostToDevice));
		
		mute_directwave<<<(ng+511)/512, 512>>>(p_cal, &d_sz[is], &d_sx[is], &d_sy[is], d_gz, d_gx, d_gy, ng, nt, dt, fm, dz, dx, dy, vmute);
		
		CHECK(cudaMemcpy(p_IO, p_cal, nt*ng*sizeof(float), cudaMemcpyDeviceToHost));
		sprintf(filename, "./data/Output/Shot/%dshot_mute.dat", is+1);
		write_data(filename, p_IO, nt*ng);
	/*************** initialize arrays with 0 ***********************************/
		CHECK(cudaMemset(g_P0,	0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(g_P1,	0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(g_P2,	0,		nnz*nnx*nny*sizeof(float)));
	/****************************** Backward ************************************/
		sprintf(filename, "./data/Output/bndr.dat");
		fp_bndr = fopen(filename, "rb");
		for(it=nt-1; it>-1; it--)
		{
		/*********************** source wavefield *******************************/
			ptr=s_P0;	s_P0=s_P1;	s_P1=s_P2;	s_P2=ptr;
			
			long offset = (long)it*nb*sizeof(float);
			fseek(fp_bndr, offset, SEEK_SET);
			fread(h_bndr, sizeof(float), nb, fp_bndr);
			CHECK(cudaMemcpy(d_bndr, h_bndr, nb * sizeof(float), cudaMemcpyHostToDevice));
			save_bndr<<<(nb+511)/512, 512>>>(nnx, nny, nnz, nx, ny, nz, nt, s_P1, d_bndr, false);
			
			step_forward3d<<<dimg,dimb>>>(s_P0, s_P1, s_P2, d_vv, damp, dz, dx, dy, dt, nnz, nnx, nny);
			cuda_add_source_lerp<<<1,1>>>(false, s_P1, &d_wlt[it], &d_sz[is], &d_sx[is], &d_sy[is], nnz, nnx, nny);
		/********************** receivers wavefield *****************************/
			cuda_record_lerp<<<(ng+511)/512, 512>>>(g_P1, p_cal, d_gz, d_gx, d_gy, ng, it, nt, nnz, nnx, nny, false);
			step_forward3d<<<dimg,dimb>>>(g_P0, g_P1, g_P2, d_vv, damp, dz, dx, dy, dt, nnz, nnx, nny);
			ptr=g_P0;	g_P0=g_P1;	g_P1=g_P2;	g_P2=ptr;
			
			if(is==0&&(it==1000))
			{
				CHECK(cudaMemcpy(temp, s_P1, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost));
				sprintf(filename,"./data/Output/snap0_b%d.dat",it);
				write_data(filename, temp, nnz*nnx*nny);
			}	
		/*********************** imaging condition ******************************/
		//	cuda_cal_illum<<<(nz*nx*ny+511)/512,512>>>(illum, g_P1, nz, nx, ny);
			cross_correlate<<<(nz*nx*ny+511)/512,512>>>(mig, s_P1, g_P1, nx, ny, nz);
		}//it loop end
		fclose(fp_bndr);
	/********************************  image  ***************************************/	
		cudaMemcpy(tmp, illum, nz*nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
		sprintf(filename,"./data/Output/illum/illum%d.dat", is+1);
		write_data(filename, tmp, nz*nx*ny);
		
		cudaMemcpy(tmp, mig, nz*nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
		sprintf(filename,"./data/Output/Mig/mig%d.dat", is+1);
		write_data(filename, tmp, nz*nx*ny);
		
		cudaEventRecord(is_t1);
		cudaEventSynchronize(is_t1);
		cudaEventElapsedTime(&seconds, is_t0, is_t1);
		is_elapsed = seconds/1000.0f;
		printPctBar(is, ns, is_elapsed);
		usleep(50000);// Sleep for 50 milliseconds
	}//is loop end
	cudaEventRecord(ns_t1);
	cudaEventSynchronize(ns_t1);
	cudaEventElapsedTime(&seconds, ns_t0, ns_t1);
	printf("\n");
	printf("Cal mig, total %d shots: %f (s)\n", ns, seconds/1000.0f);
/********************************  image  ***************************************/
	cuda_imaging<<<(nz*nx*ny+511)/512,512>>>(mig, illum, nx, ny, nz);
/***************************** filter and output file  **************************/
	sprintf(filename,"./data/Output/illum.dat");
	fp=fopen(filename,"wb");
	CHECK(cudaMemcpy(tmp, illum, nz*nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	fwrite(&tmp[0], sizeof(float), nz*nx*ny, fp);
	fclose(fp);
		
	sprintf(filename,"./data/Output/mig.dat");
	fp=fopen(filename,"wb");
	CHECK(cudaMemcpy(tmp, mig, nz*nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	fwrite(&tmp[0], sizeof(float), nz*nx*ny, fp);
	fclose(fp);
/***********************************************************************************************/
	free(tmp);		free(temp);
	free(p_IO);		free(h_bndr);
	free(gx_h);		free(gy_h);		free(gz_h);
	
	CHECK(cudaFree(p_cal));		CHECK(cudaFree(d_bndr));
	CHECK(cudaFree(d_gx));		CHECK(cudaFree(d_gy));		CHECK(cudaFree(d_gz));
	CHECK(cudaFree(s_P0));		CHECK(cudaFree(s_P1));		CHECK(cudaFree(s_P2));
	CHECK(cudaFree(g_P0));		CHECK(cudaFree(g_P1));		CHECK(cudaFree(g_P2));
	CHECK(cudaFree(illum));		CHECK(cudaFree(mig));
	
}

void rtm_ck(float *d_vv, float *damp, float *d_wlt, float *d_sz, float *d_sx, float *d_sy, char *Recefile, 
		float dz, float dx, float dy, float dt, float fm, float vmute, int nt, int ns, int ng, 
		int nz, int nx, int ny, int nnz, int nnx, int nny)
{
	int it, is;
	float *gz_h, *gx_h, *gy_h;
	float *p_IO, *tmp, *temp, *ptr;
	float *p_cal;
	float *d_gz, *d_gx, *d_gy;
	float *P0, *P1, *P2, *s_P;
	float *illum, *mig;
	
	
	cudaEvent_t is_t0, is_t1, ns_t0, ns_t1;
	cudaEventCreate(&is_t0);	cudaEventCreate(&is_t1);
	cudaEventCreate(&ns_t0);	cudaEventCreate(&ns_t1);
	float seconds = 0;
	float is_elapsed;
	
	char filename[250], wavefile[250];
	FILE *fp;
/***********************************************************************************************/
	int nk = 10;
/********************************** allocate CPU arrays ****************************************/
	tmp					= (float*)malloc(nz*nx*ny*sizeof(float));
	temp				= (float*)malloc(nnz*nnx*nny*sizeof(float));
	p_IO				= (float*)malloc(nt*ng*sizeof(float));
	
	gz_h				= (float*)malloc(ng*sizeof(float));
	gx_h				= (float*)malloc(ng*sizeof(float));
	gy_h				= (float*)malloc(ng*sizeof(float));
/******************************* allocate GPU arrays *******************************************/
	cudaSetDevice(0);// initialize device, default device=0;
	check_gpu_error("Failed to initialize device!");
	
	dim3 dimb(8, 8, 8);
	dim3 dimg((nnz+7)/8, (nnx+7)/8, (nny+7)/8);
/******************************* allocate GPU arrays *******************************************/
//*****< forward & backward & receivers wavefield >
	CHECK(cudaMalloc(&P0,			nnz*nnx*nny*sizeof(float)));
	CHECK(cudaMalloc(&P1,			nnz*nnx*nny*sizeof(float)));
	CHECK(cudaMalloc(&P2,			nnz*nnx*nny*sizeof(float)));
	CHECK(cudaMalloc(&s_P,			nz*nx*ny*sizeof(float)));
	
	CHECK(cudaMalloc(&p_cal,		nt*ng*sizeof(float)));
	
	CHECK(cudaMalloc(&d_gz,			ng*sizeof(float)));
	CHECK(cudaMalloc(&d_gx,			ng*sizeof(float)));
	CHECK(cudaMalloc(&d_gy,			ng*sizeof(float)));
//*****< The is & ns gradient ,lighting matrix >
	CHECK(cudaMalloc(&illum,		nz*nx*ny*sizeof(float)));
	CHECK(cudaMalloc(&mig,			nz*nx*ny*sizeof(float)));
/***********************************************************************************************/
	CHECK(cudaMemset(illum,		0,		nz*nx*ny*sizeof(float)));
	CHECK(cudaMemset(mig,		0,		nz*nx*ny*sizeof(float)));
	printf("===   Reverse-time Migration is start\n");
	cudaEventRecord(ns_t0);
/**************************** begin ns loop  *************************************/
	for(is=0; is<300; is++)
	{
		cudaEventRecord(is_t0);
		if(is>=0)
		{
			for(it=0; it<nt; it+=nk)
			{
				sprintf(wavefile, "./data/Output/Wavefield/check_wavefield_%d.dat", it);
				remove(wavefile);
			}
		}
		CHECK(cudaMemset(P0,	0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(P1,	0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(P2,	0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(p_cal,	0,		nt*ng*sizeof(float)));
	/*********************  set  receivers point positions ***************************/
		read_receivers(Recefile, gz_h, gx_h, gy_h, ns, is, ng, dz, dx, dy);
		CHECK(cudaMemcpy(d_gz, gz_h, ng*sizeof(float), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_gx, gx_h, ng*sizeof(float), cudaMemcpyHostToDevice));
		CHECK(cudaMemcpy(d_gy, gy_h, ng*sizeof(float), cudaMemcpyHostToDevice));
	/****************************** Forward *************************************/
		for(it=0; it<nt; it++)
		{
			cuda_add_source_lerp<<<1,1>>>(true, P1, &d_wlt[it], &d_sz[is], &d_sx[is], &d_sy[is], nnz, nnx, nny);
			step_forward3d<<<dimg,dimb>>>(P0, P1, P2, d_vv, damp, dz, dx, dy, dt, nnz, nnx, nny);
			ptr=P0;	P0=P1;	P1=P2;	P2=ptr;
		/**************** write boundary wavefield to arrays ********************/
			save_wavefield(wavefile, temp, tmp, P1, it, nt, nk, nz, nx, ny, nnz, nnx, nny);
		/************************** illumination ********************************/
			cuda_cal_illum<<<(nz*nx*ny+511)/512,512>>>(illum, P1, nz, nx, ny);
			cuda_record_lerp<<<(ng+511)/512, 512>>>(P1, p_cal, d_gz, d_gx, d_gy, ng, it, nt, nnz, nnx, nny, true);
			
			if(is==0&&it==1000)
			{
				sprintf(filename,"./data/Output/snapf.dat");
				fp=fopen(filename,"wb");
				cudaMemcpy(temp, P1, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
				fwrite(&temp[0], sizeof(float), nnz*nnx*nny, fp);
				fclose(fp);
			}		
		}//it loop end	
		sprintf(filename, "./data/Input/Shot/%dshot_obs.dat", is+1);
		read_data(filename, p_IO, nt*ng);
		CHECK(cudaMemcpy(p_cal, p_IO, nt*ng*sizeof(float), cudaMemcpyHostToDevice));
		
		mute_directwave<<<(ng+511)/512, 512>>>(p_cal, &d_sz[is], &d_sx[is], &d_sy[is], d_gz, d_gx, d_gy, ng, nt, dt, fm, dz, dx, dy, vmute);
		if(is==0 || is==ns)
		{
		CHECK(cudaMemcpy(p_IO, p_cal, nt*ng*sizeof(float), cudaMemcpyDeviceToHost));
		sprintf(filename, "./data/Output/Shot/%dshot_mute.dat", is+1);
		write_data(filename, p_IO, nt*ng);
		}
	/*************** initialize arrays with 0 ***********************************/
		CHECK(cudaMemset(P0,	0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(P1,	0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(P2,	0,		nnz*nnx*nny*sizeof(float)));
		CHECK(cudaMemset(s_P,	0,		nz*nx*ny*sizeof(float)));
	/****************************** Backward ************************************/
		for(it=nt-1; it>-1; it--)
		{
		/********************** receivers wavefield *****************************/
			cuda_record_lerp<<<(ng+511)/512, 512>>>(P1, p_cal, d_gz, d_gx, d_gy, ng, it, nt, nnz, nnx, nny, false);
			step_forward3d<<<dimg,dimb>>>(P0, P1, P2, d_vv, damp, dz, dx, dy, dt, nnz, nnx, nny);
			ptr=P0;	P0=P1;	P1=P2;	P2=ptr;
		/**************** read boundary wavefield to arrays ********************/
			load_wavefield(wavefile, tmp, P1, s_P, mig, it, nt, nk, nz, nx, ny, nnz, nnx, nny);
			
			if(is==0&&(it==1000))
			{
				CHECK(cudaMemcpy(temp, P1, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost));
				sprintf(filename,"./data/Output/snap0_b%d.dat",it);
				write_data(filename, temp, nnz*nnx*nny);
			}
		/*********************** imaging condition ******************************/
		//	cuda_cal_illum<<<(nz*nx*ny+511)/512,512>>>(illum, g_P1, nz, nx, ny);
		//	cross_correlate<<<(nz*nx*ny+511)/512,512>>>(mig, s_P1, g_P1, nx, ny, nz);
		}//it loop end
	/********************************  image  ***************************************/	
		cudaMemcpy(tmp, illum, nz*nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
		sprintf(filename,"./data/Output/illum/illum%d.dat", is+1);
		write_data(filename, tmp, nz*nx*ny);
		
		cudaMemcpy(tmp, mig, nz*nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
		sprintf(filename,"./data/Output/Mig/mig%d.dat", is+1);
		write_data(filename, tmp, nz*nx*ny);
		
		cudaEventRecord(is_t1);
		cudaEventSynchronize(is_t1);
		cudaEventElapsedTime(&seconds, is_t0, is_t1);
		is_elapsed = seconds/1000.0f;
		printPctBar(is, ns, is_elapsed);
		usleep(50000);// Sleep for 50 milliseconds
		
		remove(wavefile);
	}//is loop end
	cudaEventRecord(ns_t1);
	cudaEventSynchronize(ns_t1);
	cudaEventElapsedTime(&seconds, ns_t0, ns_t1);
	printf("\n");
	printf("Cal mig, total %d shots: %f (s)\n", ns, seconds/1000.0f);
/********************************  image  ***************************************/
	cuda_imaging<<<(nz*nx*ny+511)/512,512>>>(mig, illum, nx, ny, nz);
/***************************** filter and output file  **************************/
	sprintf(filename,"./data/Output/illum.dat");
	fp=fopen(filename,"wb");
	CHECK(cudaMemcpy(tmp, illum, nz*nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	fwrite(&tmp[0], sizeof(float), nz*nx*ny, fp);
	fclose(fp);
		
	sprintf(filename,"./data/Output/mig.dat");
	fp=fopen(filename,"wb");
	CHECK(cudaMemcpy(tmp, mig, nz*nx*ny*sizeof(float), cudaMemcpyDeviceToHost));
	fwrite(&tmp[0], sizeof(float), nz*nx*ny, fp);
	fclose(fp);
/***********************************************************************************************/
	free(tmp);		free(temp);
	free(p_IO);		
	free(gx_h);		free(gy_h);		free(gz_h);
	
	CHECK(cudaFree(p_cal));
	CHECK(cudaFree(d_gx));		CHECK(cudaFree(d_gy));		CHECK(cudaFree(d_gz));
	CHECK(cudaFree(P0));		CHECK(cudaFree(P1));		CHECK(cudaFree(P2));
	CHECK(cudaFree(illum));		CHECK(cudaFree(mig));
	
}


