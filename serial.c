#define _GNU_SOURCE

#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include <memory.h>
#include<math.h>
#include<float.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))
#define PI_POWERED M_PI*M_PI

void axmb(const int m, const int n, double uu[], double ww[], const double xg[], const double yg[]){
	int i,j;
	const double hh = xg[1] - xg[0], hh2 = hh*hh;
	memset(ww, 0, sizeof(ww));

	for(i=1;i<m-1;i++){
		for(j=0;j<n;j++){
			ww[j+n*i] += (-2.L*uu[j+n*i]+uu[j+n*(i-1)] + uu[j+n*(i+1)])/hh2;
		}
	}
	const double hh_2 = yg[1] - yg[0], hh2_2 = hh_2*hh_2;
	for(i=1;i<m-1;i++){
		for(j=1;j<n-1;j++){
			ww[j+n*i] += (-2.L*uu[j+n*i]+uu[j-1+n*i] + uu[j+1+n*i])/hh2_2;
		}
	}
	double xcos[m], ycos[n];
	for (i = 0; i < m; ++i) {
		xcos[i] = cos(M_PI*xg[i]);
	}
	for (i = 0; i < n; ++i) {
		ycos[i] = cos(M_PI*yg[i]);
	}
	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			ww[j+n*i] -= (-2.l*PI_POWERED*xcos[i]*ycos[j]);
		}
	}
}
inline void vndb(const int m, const int n, double uu[]){
	int i,j;
	for(j=0;j<n;j++) uu[j] = uu[j+1*n];
	for(j=0;j<n;j++) uu[j+n*(m-1)] = uu[j+n*(m-2)];
	for(i=0;i<m;i++) uu[n*i] = uu[1+n*i];
	for(i=0;i<m;i++) uu[n-1+n*i] = uu[n-2+n*i];
}

#define m 1000
#define n 1000

int main(int argc, char **argv){
	int i,j,k,miter,iter;
	double test0;

	double xg[m];
	double yg[n];
	double uu[m*n];
	double vv[m*n];
	double ww[m*n];

	for(i=0;i<m;i++) xg[i] = (1.L)/(double)(m-1)*(double)(i);
	for(j=0;j<n;j++) yg[j] = (1.L)/(double)(n-1)*(double)(j);

	const double hh = xg[1] - xg[0];
	const double hh2 = hh * hh;
	memset(uu, 0, sizeof(uu));
	vndb(m,n,uu);
	memcpy(vv, uu, sizeof(vv));
	miter = 10;

	double ww_[m*n];
	for (unsigned int i = 0; i < m*n; ++i) {
		ww_[i] = ww_[i]*(hh2/4.L)*0.9L;
	}

	for(iter=1;iter<=miter;iter++){
		vndb(m,n,uu);
		axmb(m,n,uu,ww,xg,yg);
		vndb(m,n,ww);
		for(i=1;i<m-1;i++){
			for(j=1;j<n-1;j++){
				uu[j+n*i] = vv[j+n*i] + ww[j+n*i]*(hh2/4.L)*0.9L;
			}
		}
		test0 = -DBL_MAX;
		double amax = -DBL_MAX;
		for (i=0;i<m*n;++i) {
			amax = max(fabs(uu[i])+1.E-8L,amax);
		}
		for(i=0;i<m*n;i++){
			test0 = max(fabs(vv[i]-uu[i]),test0);
		}
		test0 = test0/amax;
		memcpy(vv, uu, sizeof(vv));
		
		if((iter%1000) == 1 || iter < 1000) printf("%d %24.15e\n",iter, test0);
		if(test0<1.E-6L ) break;
	}

	printf("%24.15g %24.15g\n",uu[0],uu[m-1]);
	printf("%24.15g %24.15g\n",uu[n-1],uu[n-1+n*(m-1)]);
	FILE *const wp = fopen("fort.11", "w");
	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			fprintf(wp,"%g %g %g\n",xg[i],yg[i],uu[j+n*i]);
		}
		fprintf(wp,"\n");
	}
	fprintf(wp,"\n");
}
