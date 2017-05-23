
//////HW12.cpp		 
//////Author:	YUNG-SHAN SU 
//////Date: 5/22


#include<stdio.h>
#include<stdlib.h>
#include<string.h> 
#include"MAT.h"
#include<math.h>
#include<sys/time.h>
#define R 1
#define L 1
#define C 1
#define V 1

VEC function2(VEC x_t_plus_h, VEC x_t,double step);
MAT Jacobi2(VEC x_t_plus_h, VEC x_t,double step);
int main(int argc,char **argv){
	VEC x_initial(3);
	x_initial[0]=1;	x_initial[1]=0;	x_initial[2]=0;
	double tstep;
	double start=0;
	double end=10;

	
	
	///////////////Question 1.1
	printf("forwardEuler \n");
	tstep=0.1;
	MAT forwardEulerMat1 = forwardEuler(x_initial,start, end,  tstep);
	//printf("V1\tV2\tI\n");
	//forwardEulerMat1.print();
	
	VEC forwardmax1 =	columnMax(forwardEulerMat1);
	VEC forwardmin1 =	columnMin(forwardEulerMat1);
	printf("V1 Max\tV2 Max\tI_Max\n");
	forwardmax1.print();
	printf("V1 Min\tV2 Min\tI_Min\n");
	forwardmin1.print();
	printf("\n");
	
	///////////////Question 1.2
	tstep=0.01;
	MAT forwardEulerMat2 = forwardEuler(x_initial,start, end,  tstep);
	//printf("V1\tV2\tI\n");
	//forwardEulerMat2.print();
	
	
	VEC forwardmax2 =	columnMax(forwardEulerMat2);
	VEC forwardmin2 =	columnMin(forwardEulerMat2);
	printf("V1 Max\tV2 Max\tI_Max\n");
	forwardmax2.print();
	printf("V1 Min\tV2 Min\tI_Min\n");
	forwardmin2.print();
	printf("\n");
	
	///////////////Question 2.1
	printf("BackwardEuler \n");
	tstep=0.1;
	MAT backwardEulerMat1 = backwardEuler(x_initial,start, end,  tstep);
	//printf("V1\tV2\tI\n");
	//backwardEulerMat1.print();
	
	VEC backwardmax1 =	columnMax(backwardEulerMat1);
	VEC backwardmin1 =	columnMin(backwardEulerMat1);
	printf("V1 Max\tV2 Max\tI_Max\n");
	backwardmax1.print();
	printf("V1 Min\tV2 Min\tI_Min\n");
	backwardmin1.print();
	printf("\n");
	///////////////Question 2.2
	tstep=0.01;
	MAT backwardEulerMat2 = backwardEuler(x_initial,start, end,  tstep);
	//printf("V1\tV2\tI\n");
	//backwardEulerMat2.print();
	
	VEC backwardmax2 =	columnMax(backwardEulerMat2);
	VEC backwardmin2 =	columnMin(backwardEulerMat2);
	printf("V1 Max\tV2 Max\tI_Max\n");
	backwardmax2.print();
	printf("V1 Min\tV2 Min\tI_Min\n");
	backwardmin2.print();
	printf("\n");
		///////////////Question 3.1
	printf("Trep\n");
	tstep=0.1;
	MAT trepMat1 = trep(x_initial,start, end,  tstep);
	//printf("V1\tV2\tI\n");
	//trepMat1.print();
	
	VEC trepmax1 =	columnMax(trepMat1);
	VEC trepmin1 =	columnMin(trepMat1);
	printf("V1 Max\tV2 Max\tI_Max\n");
	trepmax1.print();
	printf("V1 Min\tV2 Min\tI_Min\n");
	trepmin1.print();
	printf("\n");
	///////////////Question 3.2
	tstep=0.01;
	MAT trepMat2 = trep(x_initial,start, end,  tstep);
	//printf("V1\tV2\tI\n");
	//trepMat2.print();
	
	VEC trepmax2 =	columnMax(trepMat2);
	VEC trepmin2 =	columnMin(trepMat2);
	printf("V1 Max\tV2 Max\tI_Max\n");
	trepmax2.print();
	printf("V1 Min\tV2 Min\tI_Min\n");
	trepmin2.print();
	printf("\n");
	return 0;
}	
VEC function12_forwardEuler(VEC &x_t,double step){
	VEC x_t_plus_h(x_t.len());
	
	x_t_plus_h[1]= x_t[1] + step * (x_t[2]/C);
	x_t_plus_h[2]= x_t[2] + step * ((x_t[0]-x_t[1])/L);
	x_t_plus_h[0]= -R * x_t_plus_h[2] + V; 
	return x_t_plus_h;
}
VEC function12_backwardEuler(VEC &x_t,double step){
	VEC x_t_plus_h(x_t.len());
	x_t_plus_h[2]=(x_t[2] + step*V/L - step*x_t[1]/L) / (1 + step*R/L + step*step/L/C );
	x_t_plus_h[1]=x_t[1] + step * x_t_plus_h[2]/C;
	x_t_plus_h[0]=-x_t_plus_h[2]*R+V; 
	return x_t_plus_h;
}

VEC function12_trep(VEC &x_t,double step){
	VEC x_t_plus_h(x_t.len());
	
	x_t_plus_h[2]=(x_t[2] + step*V/2/L - 2*step*x_t[1]/2/L - step*step*x_t[2]/4/L/C +x_t[0]*step/2/L) / (1 + step*R/2/L + step*step/4/L/C); 
	x_t_plus_h[1]=x_t[1] + step*(x_t_plus_h[2]+x_t[2])/2/C;
	x_t_plus_h[0]=-x_t_plus_h[2]*R+V; 
	return x_t_plus_h;
}

