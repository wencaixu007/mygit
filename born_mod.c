#include "born_mod.h"
/******************************/
void born_modeling(char FN1[],float **vp,float **ref,float **rho,int nx,int nz,
	int npd,int nt,int ndtt,int mm,int xsn,int zsn,float dx,float dz,
	float dt,float d0,float fmain,float pfac,int zrec,float **cal)
{
	int it,nxx,nzz,i,j;
	float *coeff_x1,*coeff_z1,*coeff_x2,*coeff_z2,*s;
	float **u,**w,**p,**px,**pz,**c;
	float **u0,**w0,**p0,**px0,**pz0;
	FILE *fp;
	//fp=fopen(FN1,"wb");	
	nxx=nx+2*npd;nzz=nz+2*npd;
	/

		//if((it+1)%500==0)printf("it=%d\n",it+1);
		update_p(u0,w0,p0,px0,pz0,vp,c,nxx,nzz,mm,dx,dz,dt,
	             coeff_x1,coeff_x2,coeff_z1,coeff_z2,1);
	    p0[xsn][zsn]=p0[xsn][zsn]+s[it]*vp[xsn][zsn]*dt;
	    update_vel(u0,w0,p0,rho,c,nxx,nzz,mm,dx,dz,dt,
	             coeff_x1,coeff_x2,coeff_z1,coeff_z2,1);
	             
		update_p(u,w,p,px,pz,vp,c,nxx,nzz,mm,dx,dz,dt,
	             coeff_x1,coeff_x2,coeff_z1,coeff_z2,1);
	    for(i=1;i<nxx;i++)
	    	for(j=1;j<nzz;j++)
	    	    p[i][j]=p[i][j]+ref[i][j]*p0[i][j]*vp[i][j]*dt;
	    update_vel(u,w,p,rho,c,nxx,nzz,mm,dx,dz,dt,
	             coeff_x1,coeff_x2,coeff_z1,coeff_z2,1);	             
	    /*if((it+1)%100==0) 
	    {
	    	fseek(fp,((it+1)/100-1)*nx*nz*4L,0);
	    	for(i=0;i<nx;i++)
	    		for(j=0;j<nz;j++)
	    			fwrite(&p[i+npd][j+npd],sizeof(float),1,fp);	    	
	    }*/
	    if(it%ndtt==0)
	    {
	    	for(i=0;i<nx;i++)
	    		cal[i][it/ndtt]=p[i+npd][zrec+npd-1];
	    }    	
	}
	//fclose(fp);
	free1d(coeff_x1);free1d(coeff_x2);free1d(coeff_z1);free1d(coeff_z2);
	free2d(u);free2d(w);free2d(p);free2d(px);free2d(pz);free2d(c);
	free2d(u0);free2d(w0);free2d(p0);free2d(px0);free2d(pz0);free1d(s);			
}
//################
