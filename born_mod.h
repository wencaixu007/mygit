#ifndef _BORN
#define _BORN
//@@@@@@@@@@@@@@@@@@@@@@@@@@@
#include <stdio.h>
#include "alloc.h"
#include "fd_coeff.h"
#include "io.h"
#include "pml.h"
#include "source.h"
#include "update.h"
void born_modeling(char FN1[],float **vp,float **ref,float **rho,int nx,int nz,
	int npd,int nt,int ndtt,int mm,int xsn,int zsn,float dx,float dz,
	float dt,float d0,float fmain,float pfac,int zrec,float **cal);
//@@@@@@@@@@@@@@@@@@@@@@@@@@@
#endif

