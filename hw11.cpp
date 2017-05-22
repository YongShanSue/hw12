
//////HW10.cpp		 
//////Author:	YUNG-SHAN SU 
//////Date: 5/13 




#include<stdio.h>
#include<stdlib.h>
#include<string.h> 
#include"MAT.h"
#include<math.h>
#include<sys/time.h>
void NewtonMethod(VEC& x, double errorMax,double vdd, int problem);
VEC function1(VEC x,double vdd);
MAT Jacobi1(VEC x,double vdd);
VEC function2(VEC x,double vdd);
MAT Jacobi2(VEC x,double vdd);
int main(int argc,char **argv){
	double errorMax=0.000000001;
	int dim=2;
	double vdd=1;
	double Ir;
	double I_D1;
	double I_D2;
	double I_D3;
	double I_D4;
	double r2=0.01;
	double Is=1;
	double fi=0.026;
	double k=2;
	VEC x(2);

	
	////////////////////Question1////////////////////////
	
	printf("vdd\tv1\tv2\tIr\tI_D1\tI_D2\tI_D3\tI_D4\n");
	for(double i=-100;i<=100; i++){	
		vdd =i*0.01;
		
		x[0]=vdd/2;
		x[1]=vdd/4;
	
		NewtonMethod( x, errorMax,vdd,1);
		printf("%g\t",vdd);
		x.print();
		Ir=(x[0]-x[1])/r2;
		I_D1=Is*( exp((vdd-x[0])/fi) -1   );
		I_D2= Is*(exp((-x[0])/fi)-1);
		I_D3= Is*(exp((x[1]-vdd)/fi)-1);
		I_D4= Is*( exp((x[1])/fi) -1   );
		printf("%g\t%g\t%g\t%g\t%g\n",Ir,I_D1,I_D2,I_D3,I_D4);
		
	}
	
	/*
	////////////////////Question2////////////////////////
	printf("vdd\tv1\tv2\tT1\tT2\tT3\tT4\tIr\tI_D1\tI_D2\tI_D3\tI_D4\n");
	
	VEC x_T(6);
	x_T[0]=0;
	x_T[1]=0;
	x_T[2]=300;
	x_T[3]=300;
	x_T[4]=300;
	x_T[5]=300;
	for(double i=0;i<=100; i++){	
		vdd =i*0.01;
		

		printf("%g\t",vdd);
		//printf("X\n");
		//x_T.print();
		NewtonMethod( x_T, errorMax,vdd,2,10000);
		double fi1=fi*x_T[2] /300;
		double fi2=fi*x_T[3] /300;
		double fi3=fi*x_T[4] /300;
		double fi4=fi*x_T[5] /300;
		
		
		x_T.print();
		Ir=(x_T[0]-x_T[1])/r2;
		I_D1=Is*( exp((vdd-x_T[0])/fi1) -1   );
		I_D2= Is*(exp((-x_T[0])/fi2)-1);
		I_D3= Is*(exp((x_T[1]-vdd)/fi3)-1);
		I_D4= Is*( exp((x_T[1])/fi4) -1   );
		printf("%g\t%g\t%g\t%g\t%g\n",Ir,I_D1,I_D2,I_D3,I_D4);
		
	}
	printf("vdd\tv1\tv2\tT1\tT2\tT3\tT4\tIr\tI_D1\tI_D2\tI_D3\tI_D4\n");
	x_T[0]=0;
	x_T[1]=0;
	x_T[2]=300;
	x_T[3]=300;
	x_T[4]=300;
	x_T[5]=300;
	for(double i=0;i>=-100; i--){	
		vdd =i*0.01;
		

		printf("%g\t",vdd);
		//printf("X\n");
		//x_T.print();
		NewtonMethod( x_T, errorMax,vdd,2,10000);
		double fi1=fi*x_T[2] /300;
		double fi2=fi*x_T[3] /300;
		double fi3=fi*x_T[4] /300;
		double fi4=fi*x_T[5] /300;
		
		
		x_T.print();
		Ir=(x_T[0]-x_T[1])/r2;
		I_D1=Is*( exp((vdd-x_T[0])/fi1) -1   );
		I_D2= Is*(exp((-x_T[0])/fi2)-1);
		I_D3= Is*(exp((x_T[1]-vdd)/fi3)-1);
		I_D4= Is*( exp((x_T[1])/fi4) -1   );
		printf("%g\t%g\t%g\t%g\t%g\n",Ir,I_D1,I_D2,I_D3,I_D4);
		
	}
	*/
	return 0;
}	
VEC function1(VEC x,double vdd){
	double r2=0.01;
	double Is=1;
	double fi=0.026;
	VEC Fx(x.len());
	Fx[0]=(x[0]-x[1])/r2 - Is*(exp((vdd-x[0])/fi)-1.0) - Is*(exp((-x[0])/fi)-1.0);
	Fx[1]=(x[0]-x[1])/r2 - Is*(exp(x[1]/fi)-1.0) - Is*(exp((x[1]-vdd)/fi)-1.0);

	return Fx;
}
MAT Jacobi1(VEC x,double vdd){
	double r2=0.01;
	double Is=1;
	double fi=0.026;
	double h=0.00001;
	double Fx1;
	double Fx2;
	MAT Jf(x.len());
	Fx1=(x[0]-x[1])/r2-Is*(exp((vdd-x[0])/fi)-1.0)- Is*(exp((-x[0])/fi)-1.0);
	Fx2=((x[0]+h)-x[1])/r2-Is*(exp((vdd-(x[0]+h))/fi)-1.0)- Is*(exp((-(x[0]+h))/fi)-1.0);
	Jf[0][0]=(Fx2-Fx1)/h;

	Fx1=(x[0]-x[1])/r2-Is*(exp((vdd-x[0])/fi)-1.0)- Is*(exp((-x[0])/fi)-1.0);
	Fx2=(x[0]-(x[1]+h))/r2-Is*(exp((vdd-x[0])/fi)-1)- Is*(exp((-x[0])/fi)-1.0);
	Jf[0][1]=(Fx2-Fx1)/h;
	
	Fx1=(x[0]-x[1])/r2 - Is*(exp(x[1]/fi)-1) - Is*(exp((x[1]-vdd)/fi)-1);
	Fx2=((x[0]+h)-x[1])/r2-Is*(exp(x[1]/fi)-1) - Is*(exp((x[1]-vdd)/fi)-1);
	Jf[1][0]=(Fx2-Fx1)/h;

	Fx1=(x[0]-x[1])/r2-Is*(exp(x[1]/fi)-1) - Is*(exp((x[1]-vdd)/fi)-1);
	Fx2=(x[0]-(x[1]+h))/r2-Is*(exp((x[1]+h)/fi)-1) - Is*(exp(((x[1]+h)-vdd)/fi)-1);
	Jf[1][1]=(Fx2-Fx1)/h;


	return Jf;
}
VEC function2(VEC x,double vdd){
	double r2=0.01;
	double Is=1;
	double fi0=0.026;
	double fi1;
	double fi2;
	double fi3;
	double fi4;
	double k=2;
	VEC Fx(x.len());
	double id;
	
	fi1=fi0*x[2] /300.0;
	fi2=fi0*x[3] /300.0;
	fi3=fi0*x[4] /300.0;
	fi4=fi0*x[5] /300.0;
	
	Fx[0]=(x[0]-x[1])/r2 - Is*(exp((vdd-x[0])/fi1)-1.0) - Is*(exp((-x[0])/fi2)-1.0);
	Fx[1]=(x[0]-x[1])/r2 - Is*(exp(x[1]/fi4)-1.0) - Is*(exp((x[1]-vdd)/fi3)-1.0);
	Fx[2]=x[2]-(300 + k*(vdd-x[0])* Is*(exp((vdd-x[0])/fi1)-1.0)  );
	Fx[3]=x[3]-(300 + k*(-x[0]) *  Is*(exp((-x[0])/fi2)-1.0));
	Fx[4]=x[4]-(300 + k * (x[1]-vdd) * Is*(exp((x[1]-vdd)/fi3)-1.0));
	Fx[5]=x[5]-(300 + k * (x[1]) * Is*(exp(x[1]/fi4)-1.0));
	//Fx =Fx*-1;
	return Fx;
}
MAT Jacobi2(VEC x,double vdd){
	double r2=0.01;
	double Is=1;
	double fi0=0.026;
	double h=0.00001;
	double k=2;
	double fi1;
	double fi2;
	double fi3;
	double fi4;
	VEC Fx(x.len());
	VEC Fx2(x.len());
	double fi1_diff;
	double fi2_diff;
	double fi3_diff;
	double fi4_diff;

	MAT Jf(x.len());
	fi1=fi0*x[2] /300;
	fi2=fi0*x[3] /300;
	fi3=fi0*x[4] /300;
	fi4=fi0*x[5] /300;
	
	Fx[0]=(x[0]-x[1])/r2 - Is*(exp((vdd-x[0])/fi1)-1.0) - Is*(exp((-x[0])/fi2)-1.0);
	Fx[1]=(x[0]-x[1])/r2 - Is*(exp(x[1]/fi4)-1.0) - Is*(exp((x[1]-vdd)/fi3)-1.0);
	Fx[2]=x[2]-(300 + k  * Is*(exp((vdd-x[0])/fi1)-1.0)  * (vdd-x[0]) );
	Fx[3]=x[3]-(300 + k*  Is*(exp((-x[0])/fi2)-1.0)   *(-1*x[0])  );
	Fx[4]=x[4]-(300 + k * Is*(exp((x[1]-vdd)/fi3)-1.0) * (x[1]-vdd) );
	Fx[5]=x[5]-(300 + k  * Is*( exp(x[1]/fi4) -1.0 ) * (x[1]) );
	
	
	//VEC xreg(x);
	for(int i=0;i<x.len();i++){

		for(int j=0;j<x.len();j++){
			VEC xreg(x);
			xreg[j] += h;
			fi1_diff=fi0*xreg[2] /300;
			fi2_diff=fi0*xreg[3] /300;
			fi3_diff=fi0*xreg[4] /300;
			fi4_diff=fi0*xreg[5] /300;
			if(i==0)
				Fx2[j]= (xreg[0]-xreg[1])/r2 - Is*(exp((vdd-xreg[0])/fi1_diff)-1.0) - Is*(exp((-xreg[0])/fi2_diff)-1.0);
			else if(i==1)
				Fx2[j]= (xreg[0]-xreg[1])/r2 - Is*(exp(xreg[1]/fi4_diff)-1.0) - Is*(exp((xreg[1]-vdd)/fi3_diff)-1.0);
			else if(i==2)
				Fx2[j]= xreg[2]-(300 + k*(vdd-xreg[0])* Is*(exp((vdd-xreg[0])/fi1_diff)-1.0)  );
			else if(i==3)
				Fx2[j]= xreg[3]-(300 + k*(-xreg[0]) *  Is*(exp((-xreg[0])/fi2_diff)-1.0));
			else if(i==4)
				Fx2[j]= xreg[4]-(300 + k * (xreg[1]-vdd) * Is*(exp((xreg[1]-vdd)/fi3_diff)-1.0));
			else	
				Fx2[j]= xreg[5]-(300 + k * (xreg[1]) * Is*(exp(xreg[1]/fi4_diff)-1.0));

			Jf[i][j]=(Fx2[j]-Fx[i])/h;
		}

	}
	return Jf;
}
void NewtonMethod(VEC& x, double errorMax,double vdd, int problem){
	int k=0;
	double error=errorMax+1;
	VEC Fx(x.len());
	MAT Jf(x.len());
	VEC deltaX(x.len());
	VEC reg(x.len());
	if (problem==1){
		while(error>=errorMax){

			//x.print();
			Fx=function1(x,vdd);
			//Fx.print();//printf("\n");
			Jf=Jacobi1(x,vdd);
			//Jf.print();
			getLinearSolution(Jf,deltaX,(Fx*-1));

			x = x+deltaX;
			reg=function1(x,vdd);
			error=max_norm(deltaX);
			//deltaX.print(); printf("\n");
			//printf("\n");
			//printf("%g\n",error);
		}
	}
	else{
		while(error>=errorMax){
			Fx=function2(x,vdd);

			Jf=Jacobi2(x,vdd);
			getLinearSolution(Jf,deltaX,(Fx*-1));
			x=x+deltaX;
			reg=function2(x,vdd);
			error=max_norm(reg);

		}
	}
	
			
}
