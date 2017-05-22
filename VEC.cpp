//////VEC.cpp		 
//////Author:	YUNG-SHAN SU 
//////Date: 3/19 



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "VEC.h"

VEC::VEC(int n){		// uninit constructor, val set to 0
	dim=n;
	val=(double *)calloc(n,sizeof(double));
}
VEC::VEC(const VEC &v1){	// copy constructor
	dim=v1.dim;
	val=(double*)malloc(dim*sizeof(double));
	for (int i=0;i<dim;i++)
		val[i]=v1.val[i];
}
VEC::VEC(int n,double *v){// init constructor
	dim=n;
	val=(double*)malloc(dim*sizeof(double));
	for (int i=0;i<dim;i++)
		val[i]=v[i];
}

VEC::~VEC(){// destructor
	free(val);
}

int VEC::len(){// dimension of the vector
	return(dim);
}
VEC & VEC::operator-(){ // unary operator, negative value
	for (int i=0;i<dim;i++)
		val[i]=-val[i];
	return *this;
}

VEC & VEC::operator=(const VEC &v1){ // assignment
	for (int i=0;i<v1.dim;i++)
		val[i]=v1.val[i];
	return *this;
}
VEC & VEC::operator=(const double &a){ // assignment
	for (int i=0;i<dim;i++)
		val[i]=a;
	return *this;
}
VEC & VEC::operator+=(const VEC v1){// V += v1;
	for(int i=0;i<dim;i++)
		val[i]+=v1.val[i];
	return *this;
}
VEC & VEC::operator-=(const VEC v1){ // V -= v1;
	for(int i=0;i<dim;i++)
		val[i]-=v1.val[i];
	return *this;
}
VEC & VEC::operator*=(double a){ // V *= dbl;
	for(int i=0;i<dim;i++)
		val[i]*=a;
	return *this;
}
VEC & VEC::operator/=(double a){// V /= dbl;
	for(int i=0;i<dim;i++)
		val[i]/=a;
	return *this;
}
VEC VEC::operator+(const VEC &v1){ // V + v1
	VEC v2(*this);
	for(int i=0;i<v2.dim;i++)
		v2.val[i]+=v1.val[i];
	return v2; 
}
VEC VEC::operator-(const VEC &v1){ // V - v1
	VEC v2(*this);
	for(int i=0;i<v2.dim;i++)
		v2.val[i]-=v1.val[i];
	return v2; 
}
double VEC::operator*(VEC &v1){// inner product
	double sum=0;
	for (int i=0;i<dim;i++){
        if(val[i]!=0 && v1.val[i]!=0 )
		    sum+=val[i]*v1.val[i];
            
    }
	return sum;
} 
VEC VEC::operator*(double a){// V * dbl
	VEC v2(*this);
	for(int i=0;i<dim;i++)
		v2.val[i]*=a;
	return v2;
}
VEC VEC::operator/(double a){ // V / dbl
	VEC v2(*this);
	for(int i=0;i<dim;i++)
		v2.val[i]/=a;
	if (a==0)
		printf("your division have a problem: Divisor is 0");
	return v2;
}
double & VEC::operator[](int n){ // indexing
	if (n<0) 
		n=0;
	else if (n>=dim) 
		n=dim-1;
	return val[n];
}
VEC operator*(double a,const VEC v1){
	VEC v2(v1);
	v2*=a;
	return v2;
}

VEC *newVEC(int n){ // create dynamic VEC
	VEC *vptr;
	vptr=(VEC *)malloc(sizeof(VEC));
	vptr->dim=n;
	vptr->val=(double*)calloc(n,sizeof(double));
	return vptr;
}
void	VEC::print(){
	for(int i=0;i<dim;i++)
		printf("%g\t",val[i]);
	//printf("\n");
}
VEC VEC::absolute(){
	VEC v2(*this);
	//double reg;
	for(int i=0;i<v2.len();i++)
		v2[i]=fabs(v2[i]);
	return v2;
}