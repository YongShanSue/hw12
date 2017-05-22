//////MAT.cpp		 
//////Author:	YUNG-SHAN SU 
//////Date: 4/11

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MAT.h"
#include <sys/time.h>
MAT::MAT(int dim){// uninit constructor
	n=dim;
	va=(VEC **)malloc(n*sizeof(VEC*));
	for(int i=0;i<n;i++)
		va[i]=newVEC(n);
}
MAT::MAT(const MAT &m1){ // copy constructor
	VEC **vsrc=m1.va; // to get around not indexing const MAT
	n=m1.n;
	va=(VEC **)malloc(n*sizeof(VEC*));
	for (int i=0; i<n; i++) {
		va[i]=newVEC(n);
		(*va[i])=(*vsrc[i]); // VEC assignment
	}
}
MAT::MAT(int dim,double *v){ // init constructor
	n=dim;
	va=(VEC **)malloc(n*sizeof(VEC*));
	for (int i=0; i<n; i++) {
		va[i]=newVEC(n);
		for (int j=0; j<n; j++) {
			(*va[i])[j]=*(v++); // array indexing + VEC indexing
		}
	}
}
MAT::~MAT(){ // destructor
	for (int i=n-1; i>=0; i--)
		(*va[i]).~VEC();
	free(va);
}
int MAT::dim(){ // return dimension of the matrix
	return n;
}
MAT MAT::tpose(){ // transpose
	MAT mnew(n);
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			mnew[i][j]=(*va[j])[i];
		}
	}
	return mnew;
}
MAT & MAT::operator-(){ // unary operator, negative value
	for (int i=0; i<n; i++)
		(*va[i])=-(*va[i]);
	return *this;
}
MAT &MAT::operator=(MAT m1){ // assignment
	for (int i=0; i<n; i++)
		(*va[i])=m1[i];
	return *this;
}
MAT &MAT::operator=(double a){ // assignment
	for (int i=0; i<n; i++)
		(*va[i])=a;
	return *this;
}
MAT &MAT::operator+=(MAT &m1){ // m += m1;
	for (int i=0; i<n; i++)
		(*va[i])+=m1[i];
	return *this;
}
MAT &MAT::operator-=(MAT &m1){ // m -= m1;
	for (int i=0; i<n; i++)
		(*va[i])-=m1[i];
	return *this;
}
MAT &MAT::operator*=(double a){ // m *= dbl;
	for (int i=0; i<n; i++)
		(*va[i])*=a;
	return *this;
}
MAT &MAT::operator/=(double a){ // m /= dbl;
	for (int i=0; i<n; i++)
		(*va[i])/=a;
	return *this;
}
MAT MAT::operator+(MAT m1){ // m + m1
	MAT m2(n);
	for (int i=0; i<n; i++)
		for(int j=0;j<n;j++)
			m2[i][j]=(*va[i])[j]+m1[i][j];
	return m2;
}
MAT MAT:: operator-(MAT m1){ // m - m1
	MAT m2(n);
	for (int i=0; i<n; i++)
		for(int j=0;j<n;j++)
			m2[i][j]=(*va[i])[j]-m1[i][j];
	return m2;
}
MAT MAT:: operator*(MAT m1){// m1 * m2
	MAT z(n);
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			z[i][j]=0;
			for (int k=0; k<n; k++)
				z[i][j]+=((*va[i])[k]*m1[k][j]);
		}
	}
	return z;
}
VEC &  MAT::operator[](int m){ // m'th row
	return *va[m];
}
VEC MAT::operator*(VEC &v1){ // m x v1
	VEC s(n);
	for (int i=0; i<n; i++) {
		s[i]=(*va[i])*v1; // VEC inner product
	}
	return s;
}
MAT MAT::operator*(double a){ // m * dbl
	for(int i=0;i<n;i++)
		(*va[i] ) =(*va[i])*a;
}
MAT MAT::operator/(double a){ // m / dbl
	for(int i=0;i<n;i++)
		(*va[i] ) =(*va[i])/a;
}

void MAT::print(){		//Print the matrix
	for(int row=0;row<n;row++){
		va[row]->print();
	}
		
}

MAT operator*(double a,MAT &m1){ // dbl x m
	int n=m1.dim();
	MAT z(n);
	for(int i=0;i<n;i++)
		z[i]=m1[i]*a;
	return z;
}
VEC operator*(VEC &v1,MAT &m1){ // vT x m
	VEC v2(m1.dim());
	for (int i=0; i<m1.dim(); i++) {
		v2[i]=0;
		for (int j=0; j<m1.dim(); j++) {
			v2[i] += v1[j]*m1[j][i];
		}
	}
	return v2;
}

MAT &luFact(MAT &m1){ // LU decomposition
	int i,j,k;
	for (i=0; i<m1.dim(); i++) {
		// copy a[i][j] to u[i][j] needs no action due to in-place LU
 		for (j=i+1; j<m1.dim(); j++) { // form l[j][i]
 			m1[j][i] /= m1[i][i];
 		}
 		for (j=i+1; j<m1.dim(); j++) { // update lower submatrix
 			for (k=i+1; k<m1.dim(); k++) {
 				m1[j][k] -= m1[j][i]*m1[i][k];
 			}
 		}
 	}
	return m1;
}
VEC fwdSubs(MAT &m1,VEC b){ // forward substitution
	int i,j;
	VEC y=b;		// initialize y to b
 	for (i=0; i<m1.dim(); i++)
 		for (j=i+1; j<m1.dim(); j++)
 			y[j] -= m1[j][i]*y[i];
 	return y;
}
VEC bckSubs(MAT &m1,VEC b){// backward substitution
	int i,j,k;
	VEC x=b	;		// initialize x to y
	for (i=m1.dim()-1; i>=0; i--) {
 		x[i] /= m1[i][i];
 		for (j=i-1; j>=0; j--)
 			x[j] -= m1[j][i]*x[i];
 	}
 	return x;
}
MAT  getResistorNetworksMatrix(int dim){ 	//Get the linear system network
	double resister=2000/dim;
	double conductence=1/resister;
	MAT system((dim+1)*(dim+1));
    system=0;
	VEC b((dim+1)*(dim+1));b[dim/2]=1;  //For remember
	//Process the first row//////////////////////////////////////////////////////
	
	//////Process the first column
	system[0][0]=conductence*2;
	system[0][1]=-conductence;
	system[0][0+dim+1]=-conductence;
	//////Process the second column to the last 2 column
	for(int column=1;column<dim;column++){
		system[column][column]=conductence*3;
		system[column][column-1]=-conductence;
		system[column][column+1]=-conductence;
		system[column][column+dim+1]=-conductence;
	}
	//////Process the last column
	system[dim][dim]=conductence*2;
	system[dim][dim-1]=-conductence;
	system[dim][dim+dim+1]=-conductence;
	
	//Process from the second row to the last two row//////////////////////
	for(int row=1;row<dim;row++){
		//////Process the first column
		system[row*(dim+1)][row*(dim+1)]=conductence*3;
		system[row*(dim+1)][row*(dim+1)-(dim+1)]=-conductence;
		system[row*(dim+1)][row*(dim+1)+1]=-conductence;
		system[row*(dim+1)][row*(dim+1)+(dim+1)]=-conductence;
		//////Process the second column to the last 2 column
		for(int column=1;column<dim;column++){
			system[row*(dim+1)+column][row*(dim+1)+column]=conductence*4;
			system[row*(dim+1)+column][row*(dim+1)+column-(dim+1)]=-conductence;
			system[row*(dim+1)+column][row*(dim+1)+column-1]=-conductence;
			system[row*(dim+1)+column][row*(dim+1)+column+1]=-conductence;
			system[row*(dim+1)+column][row*(dim+1)+column+(dim+1)]=-conductence;
		}
		//////Process the last column
		system[row*(dim+1)+dim][row*(dim+1)+dim]=conductence*3;
		system[row*(dim+1)+dim][row*(dim+1)+dim-(dim+1)]=-conductence;
		system[row*(dim+1)+dim][row*(dim+1)+dim-1]=-conductence;
		system[row*(dim+1)+dim][row*(dim+1)+dim+(dim+1)]=-conductence;
	}
	//Process the last row/////////////////////////////////////////
	//////Process the first column
	system[dim*(dim+1)][dim*(dim+1)]=conductence*2;
	system[dim*(dim+1)][dim*(dim+1)-(dim+1)]=-conductence;
	system[dim*(dim+1)][dim*(dim+1)+1]=-conductence;
	//////Process the second column to the last 2 column
	for(int column=1;column<dim;column++){
		system[dim*(dim+1)+column][dim*(dim+1)+column]=conductence*3;
		system[dim*(dim+1)+column][dim*(dim+1)+column-(dim+1)]=-conductence;
		system[dim*(dim+1)+column][dim*(dim+1)+column-1]=-conductence;
		system[dim*(dim+1)+column][dim*(dim+1)+column+1]=-conductence;
	}
	//////Process the last column
	system[dim*(dim+1)+dim][dim*(dim+1)+dim]=conductence*2;
	system[dim*(dim+1)+dim][dim*(dim+1)+dim-(dim+1)]=-conductence;
	system[dim*(dim+1)+dim][dim*(dim+1)+dim-1]=-conductence;
	////////////Set the Boundary condition
	system[dim/2]=0;
	system[dim/2][dim/2]=1;
	system[(dim)*(dim+1)+dim/2]=0;
	system[(dim)*(dim+1)+dim/2][(dim)*(dim+1)+dim/2]=1;
	return system;
}

//////////////////////HW 4 function///////////////////////////
int jacobi(MAT &A,VEC b,VEC &x,int maxIter,double tol){
	MAT ef=A;			//E+F
    ef*=-1;
   struct timeval start;               //Variable start time
   struct timeval stop;                //Variable stop timei
   gettimeofday(&start,NULL);
   double duration;
    MAT d_inv(A.dim());		//inverse(D)
	VEC x_reg=x;		//X_register
	VEC x_subtraction(x.len());
	double norm=0;
	d_inv=0;
    ///Set D
	for(int i=0;i<ef.dim();i++){
		ef[i][i]=0;
		d_inv[i][i]=1.0/A[i][i];
	}
	for(int i=0;i<maxIter;i++){
        x_reg=x;
        for(int j=0;j<x.len();j++)
            x[j]=d_inv[j][j]*(b[j]+(ef[j]*x_reg));
        // gettimeofday(&stop,NULL);
		x_subtraction=x-x_reg;
        //1-Norm check
		//norm=one_norm(x_subtraction);
		//2-Norm check
		//norm=two_norm(x_subtraction);
		//max-Norm check
		norm=max_norm(x_subtraction);	
        if(norm<=tol){
            gettimeofday(&stop,NULL);
            duration=(stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/1000000.0;    
            printf("CPU time:\t%g\n",duration);
			return (i+1);
        }
	}
		
         gettimeofday(&stop,NULL);
         duration=(stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/1000000.0;    
           
         printf("Duration: %g\n",duration);
	
	return maxIter;
}
int gaussSeidel(MAT &A,VEC b,VEC &x,int maxIter,double tol){
    MAT e(A.dim());
    MAT f(A.dim());
   struct timeval start;               //Variable start time
   struct timeval stop;                //Variable stop timei
   gettimeofday(&start,NULL);
   double duration;
    MAT d_inv(A.dim());
    e=0;
    f=0;
    d_inv=0;
    VEC x_reg=x;
    VEC x_subtraction(x.len());
    int circuit_dim=40;
    double norm;
    for(int row =1;row<e.dim();row++){
        for(int column=0; column<row;column++){
            e[row][column] = -1*A[row][column] ;  
        }
    }    
    for(int row=0;row<f.dim();row++){
        for(int column=row+1;column<f.dim();column++)
            f[row][column] = -1*A[row][column];
    }
    for(int i=0;i<d_inv.dim();i++)
        d_inv[i][i]=1.0/A[i][i];
    for(int i=0;i<maxIter;i++){
        for(int j=0;j<x.len();j++){
            x[j]=(b[j]+e[j]*x+f[j]*x_reg)*d_inv[j][j];
        }
        x_subtraction=x-x_reg;
        x_reg=x;
        //1-Norm check
		//norm=one_norm(x_subtraction);
		//2-Norm check
		//norm=two_norm(x_subtraction);
		//max-Norm check
		norm=max_norm(x_subtraction);	
		
		if(norm<=tol){
            gettimeofday(&stop,NULL);
            duration=(stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/1000000.0;    
            printf("CPU time:\t%g\n",duration);
			return (i+1);
        }
        
		
	}
         gettimeofday(&stop,NULL);
         duration=(stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/1000000.0;    
           
         printf("Duration: %g\n",duration);
    return maxIter;

}
int sgs(MAT &A,VEC b,VEC &x,int maxIter,double tol){
    
    MAT e(A.dim());
    MAT f(A.dim());
    MAT d_inv(A.dim());
   struct timeval start;               //Variable start time
   struct timeval stop;                //Variable stop timei
   gettimeofday(&start,NULL);
   double duration;
    e=0;
    f=0;
    d_inv=0;
    VEC x_reg=x;
    VEC x_reg_1=x;
    VEC x_subtraction(x.len());
    int circuit_dim=40;
    double norm;
    for(int row =1;row<e.dim();row++){
        for(int column=0; column<row;column++){
            e[row][column] = -1*A[row][column] ;  
        }
    }    
    for(int row=0;row<f.dim();row++){
        for(int column=row+1;column<f.dim();column++)
            f[row][column] = -1*A[row][column];
    }
    for(int i=0;i<d_inv.dim();i++)
        d_inv[i][i]=1.0/A[i][i];
    for(int i=0;i<maxIter;i++){
        for(int j=0;j<x.len();j++){
            x[j]=(b[j]+e[j]*x+f[j]*x_reg)*d_inv[j][j];
        }
        x_reg_1=x;
        for(int j=x.len()-1;j>=0;j--){
            x[j]=(b[j]+e[j]*x_reg_1+f[j]*x)*d_inv[j][j];
        }
        x_subtraction=x-x_reg;
        x_reg=x;
        //1-Norm check
		//norm=one_norm(x_subtraction);
		//2-Norm check
		//norm=two_norm(x_subtraction);
		//max-Norm check
		norm=max_norm(x_subtraction);	
		if((i+1)%100==0)
		printf("%g\t",norm);
		if(norm<=tol){
            gettimeofday(&stop,NULL);
            duration=(stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/1000000.0;    
            printf("CPU time:\t%g\n",duration);
			return (i+1);
        }
		
    }
         gettimeofday(&stop,NULL);
         duration=(stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/1000000.0;    
           
         printf("Duration: %g\n",duration);
    return maxIter;
}
/////////////////////Calculate one-norm////////////////
double one_norm(VEC &va){
    double norm=0;
    for(int j=0;j<va.len();j++)
        norm+=fabs(va[j]);
    return norm;

}
/////////////////////Calculate two-norm////////////////
double two_norm(VEC &va){
    double norm;
    norm=sqrt(va*va);
    return norm;
}
/////////////////////Calculate max-norm////////////////
double max_norm(VEC &va){
    double norm=fabs(va[0]);
    for(int j=1;j<va.len();j++){
        if(norm<fabs(va[j]))
            norm=fabs(va[j]);
    }
    return norm;

}
	
	
int cg(MAT &A,VEC b,VEC &x,int maxIter, double tol){
	struct timeval start;               //Variable start time
    struct timeval stop;                //Variable stop timei
    gettimeofday(&start,NULL);			//Record star time
    double duration;					//Variable time duration
	VEC p=b-A*x;						//Set p=b-A*x
	VEC r=p;							//Set r=p
	VEC r_reg(r.len());					//Set r(k+1)
	double a_k,b_k;						//Variable a_k, b_k
	VEC x_reg(x.len());					//Variable x(k+1)
	VEC reg(x.len());					//Varible reg=A*p
	VEC x_v_hw5(3);						//Vector x_ne,V_ea,V_sw
	int circuit_dim=(int) (sqrt(x.len())  -1) ;
	//printf("circuit_dim\t%d\n",circuit_dim);
	//printf("IterationNumber\tV_ne\tV_ea\t_V_sw\tError\n");
	for(int i=0;i<maxIter;i++){
		///Conjugate gradient method
		reg=A*p;
		a_k=(r*r)/(p*(reg));
		x_reg=x+p*a_k;
		r_reg=r-(A*p)*a_k;
		b_k=(r_reg*r_reg)/(r*r);
		p=r_reg+b_k*p;
		r=r_reg;
		x=x_reg;
		////Print the output when running 50 times
		
		if((i+1)%10==0){
			//printf("%d\t",i+1);
			x_v_hw5[0]=x[circuit_dim];
			x_v_hw5[1]=x[(circuit_dim+1)*circuit_dim/2+(circuit_dim)];
			x_v_hw5[2]=x[(circuit_dim+1)*circuit_dim];
			//x_v_hw5.print();
			//printf("%g\n",sqrt((r*r)/x.len())  );
		}
		
		////Print the output when error is low than tolerance and break the loop
		if(sqrt(  (r*r) / A.dim() )<=tol){
			//printf("%d\t",i+1);
			x_v_hw5[0]=x[circuit_dim];
			x_v_hw5[1]=x[(circuit_dim+1)*circuit_dim/2+(circuit_dim)];
			x_v_hw5[2]=x[(circuit_dim+1)*circuit_dim];
			//x_v_hw5.print();
			//printf("%g\n",sqrt((r*r)/x.len())  );
			//gettimeofday(&stop,NULL);
			//duration=(stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/1000000.0;    
   			//printf("ExecutionTime\t%g\n",duration);
			return i+1;
		}
	}
	//gettimeofday(&stop,NULL);
	//duration=(stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/1000000.0;    
    //printf("ExecutionTime\t%g\n",duration);
	return maxIter;
}
/////////HW6 function//////////////
int EVpwr(MAT &A,VEC &q0,double &lambda,double tol,int maxiter){
	VEC q1=q0;
	VEC q2(q0.len());
	VEC z(q0.len());
	VEC reg(q0.len());
	double lambda1;
	double lambda2;
	double duration;                    //Variable time duration
    struct timeval start;               //Variable start time
    struct timeval stop;                //Variable stop timei
	int count =0;
	double err=100;
	//printf("Epson1\n");
	//printf("Iteration\tError\n");
	gettimeofday(&start,NULL);			//Record star time
		////Power method algorithm
	while(tol<=err && count<maxiter){
		z=A*q1;
		q2=z/two_norm(z);
		reg=A*q2;
		lambda2=q2*(reg);
		
		err=epson_1(lambda1,lambda2);
		//err=epson_2(q1,q2);
		//err=epson_3(A,q2,lambda2);
		//err=epson_4(A,q2,lambda2);
		q1=q2;
		lambda1=lambda2;
		count++;
		//printf("%d\t%e\n",count,err);
	}
	q0=q2;
	lambda=lambda2;
	gettimeofday(&stop,NULL);			//Record star time
	duration=(stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/1000000.0;    
	printf("PWrExecutionTime\t%g\t",duration);
	//printf("%g\t",duration);
	return count;
}
int EViPwr(MAT &A,VEC &q0,double &lambda,double tol,int maxiter){
	double err=100;
	VEC z(q0.len());
	VEC reg(q0.len());
	int count =0;
	double lambda1=1;
	double lambda2=1;
	double duration;                    //Variable time duration
    struct timeval start;               //Variable start time
    struct timeval stop;                //Variable stop timei
	gettimeofday(&start,NULL);
	////Inverse Power method algorithm
	while(err>=tol && count<maxiter){
		cg(A,q0,z,300,0.0000000007);
		q0=z/two_norm(z);
		reg=A*q0;
		lambda2=q0*(reg);
		err=epson_1(lambda1,lambda2);
		lambda1=lambda2;
		count++;
	}
	lambda=lambda2;
	gettimeofday(&stop,NULL);			//Record star time
	duration=(stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/1000000.0;    
	printf("EviPWrExecutionTime\t%g\t",duration);
	//printf("%g\t",duration);
	return count;
}
int EViPwrShft(MAT &A,VEC &q0,double &lambda,double mu,double tol,int maxiter){
	MAT linearsystem=A;
	MAT In(linearsystem.dim());
	In=0;
	double lambda1=1;
	double lambda2=1;
	double duration;                    //Variable time duration
    struct timeval start;               //Variable start time
    struct timeval stop;                //Variable stop timei
    VEC z(q0.len());
	VEC reg(q0.len());
	gettimeofday(&start,NULL);
	for(int i=0;i<linearsystem.dim();i++)
		In[i][i]=1;
	linearsystem=linearsystem-(mu*In);
	double err=100;
	
	int count =0; 
	////Inverse Power method with shifting algorithm
	while(err>=tol && count<maxiter){
 		cg(linearsystem,q0,z,300,0.0000000007);
		q0=z/two_norm(z);
		reg=A*q0;
		lambda2=q0*(reg);
		err=epson_1(lambda1,lambda2);
		lambda1=lambda2;
		count++;
	}
	lambda=lambda2;
	gettimeofday(&stop,NULL);			//Record star time
	duration=(stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/1000000.0;    
	printf("EviPWrShftExecutionTime\t%g\t",duration);
	//printf("%g\t",duration);
	return count;
}

double epson_1(double &lambda1,double &lambda2){
	double err;
	err=fabs(lambda2-lambda1);
	return err;
}
double epson_2(VEC &q1,VEC &q2){
	VEC q(q1.len());
	double err;
	q=q2-q1;
	err=two_norm(q);
	return err;
}
double epson_3(MAT &A,VEC &q,double &lambda){
	VEC r(q.len());
	double err;
	r=A*q-lambda*q;
	err=two_norm(r);
	return err;

}
double epson_4(MAT &A,VEC &q,double &lambda){
	VEC r(q.len());
	VEC w(q.len());
	double err;
	r=A*q-lambda*q;
	VEC qa=q*A;
	w=(qa)/two_norm(qa);
	err=two_norm(r)/fabs(w*q);
	return err;
}


int EVqr(MAT &A,double tol,int maxiter){
	MAT r(A.dim());						//Declare r 
	MAT q(A.dim());						//Declare q
	double err=1+tol;
	r=0;
	q=0;
	int count=0; 
	while(count<maxiter&&err>tol){
		//QR iteration Algorithm
		A=A.tpose();
		r[0][0]=sqrt(A[0]*A[0]);		//r11=sqrt(A1*A1)
		q[0]=A[0]/r[0][0];				//q1=A1/r11
		for(int j=1;j<A.dim();j++){
			q[j]=A[j];
			for(int i=0;i<j;i++){
				r[j][i]=q[i]*q[j];		//rij=qi*Aj
				q[j]=q[j]-r[j][i]*q[i]; 
			}
			r[j][j]=sqrt(q[j]*q[j]);	//rjj=Aj-sigma(rij*qi)
			q[j]=(q[j])/r[j][j];		//qj=qj/rjj
		}
		A=(r.tpose())*(q.tpose());
		//Calculate error
		err=A[1][0];
		for(int i=2;i<A.dim();i++)
			if(A[i][i-1]>err)
				err=A[i][i-1];
		count++;

	}
	return count;

}
int EVqrShifted(MAT &A,double mu,double tol,int maxiter){
	MAT r(A.dim());						//Declare r 
	MAT q(A.dim());						//Declare q
	MAT In(A.dim());
	double err=1+tol;
	VEC reg(A.dim());
	In=0;
	r=0;
	q=0;
	for(int i=0;i<A.dim();i++)
		In[i][i]=1;
	int count=0; 
	//Shifted QR iteration Algorithm
	while(count<maxiter&&err>=tol){		
		A=A-(mu*In);					//A=A-mu*In
		A=A.tpose();
		r[0][0]=sqrt(A[0]*A[0]);		//r11=sqrt(A1*A1)
		q[0]=A[0]/r[0][0];				//q1=A1/r11
		for(int j=1;j<A.dim();j++){
			q[j]=A[j];
			for(int i=0;i<j;i++){
				r[i][j]=q[i]*q[j];	
				q[j]-=r[i][j]*q[i];
			}
			r[j][j]=sqrt(q[j]*q[j]);
			q[j]=(q[j])/r[j][j];		//qj=qj/rjj
		}
		A=(r)*(q.tpose())+(mu*In);//A=A+mu*In
		//Calculate error
		err=fabs(A[1][0]);
		for(int i=2;i<A.dim();i++)
			if(fabs(A[i][i-1])>err)
				err=fabs(A[i][i-1]);
		count++;
		//printf("%d\t%g\t%g\n",count,err,A[0][0]);
	}
	return count;
}
double NEV(double x,VEC &XS,VEC &YS,int n){			//Non-recursive Neville’s Algorithm
	double NS[n];
	int i,j,k;
	for (i=0; i<n; i++) 
		NS[i]=YS[i];
	for (k=1; k<n; k++) {
		for (j=0; j<n-k; j++) {
			NS[j]=((x-XS[j])*NS[j+1]-(x-XS[k+j])*NS[j])/(XS[j+k]-XS[j]);
		}
	}
	return NS[0];
}
double Lagrange(double x,VEC &XDATA,VEC &YDATA){		//Lagrange Algorithm
	int n=XDATA.len();
	return NEV( x,XDATA,YDATA,n);
}

void splineM(int N,VEC &X,VEC &Y,VEC &M){// generate spline momentum M
	MAT  system(N);
	VEC  d(N);
	system=0;
	
	/////////Boundary condition:	Usong zero boundary moments
	d[0]=0;
	d[N-1]=0;
	system[0][1]=0;		//lambda0 = 0
	system[N-1][N-2]=0; //un = 0
	/////////Set the matrix diagonal
	for(int i=0;i<N;i++)
		system[i][i]=2;
	/////////Set the matrix lambda and u and d
	for(int i=1;i<N-1;i++){
		system[i][i-1]=(X[i]-X[i-1])/(X[i+1]-X[i-1]);		//Set u
		system[i][i+1]=(X[i+1]-X[i])/(X[i+1]-X[i-1]);		//Set lambda
	}

	for(int i=1;i<N-1;i++){
		d[i]=6.0/(X[i+1]-X[i-1])*(	(Y[i+1]-Y[i])/(X[i+1]-X[i])- (Y[i]-Y[i-1]) / (X[i]-X[i-1])	);
	}
	////////////////LU decomposition////////////////////////////////////////////////
	luFact(system);			//lu decomposition
	/////////////Produce u///////////
	MAT u(N);							//Declare u
	u=0;								//Set u to 0
	for(int i=0;i<N;i++)
		for(int j=i;j<N;j++)
			u[i][j]=system[i][j];//Produce u
	/////////////Produce l///////////
	MAT l(N);							//Delare l
	l=0;								//Set l to 0
	for(int i=0;i<N;i++){
		for(int j=0;j<i;j++)
			l[i][j]=system[i][j];//Produce l
		l[i][i]=1;
	}
	VEC y=fwdSubs(l,d);
	M=bckSubs(u,y);
	return;
} 
double spline(double x,int N,VEC &X,VEC &Y,VEC &M){// spline interp at x
	int interval_index=0;
	int check=0;
	////////////Find the interpolation interval
	while(check!=1 && interval_index<=X[N-1] ){
		if (x>=X[interval_index]  && x<=X[interval_index+1])
			check=1;
		else 
			interval_index++;
	}
	double h=X[interval_index+1]-X[interval_index];
	double a=Y[interval_index];
	double b=(Y[interval_index+1]-Y[interval_index])/h-h/6*(M[interval_index+1]+2*M[interval_index]);
	double r=M[interval_index]/2;
	double delta=(M[interval_index+1]-M[interval_index])/(6*h);
	return (a + b*(x-X[interval_index]) + r* (x-X[interval_index])*(x-X[interval_index]) + delta*(x-X[interval_index])*(x-X[interval_index])*(x-X[interval_index]));
}
void getLinearSolution(MAT A,VEC & y,VEC b){
	luFact(A);			//lu decomposition
	/////////////Produce u///////////
	double N=A.dim();
	MAT u(N);							//Declare u
	u=0;								//Set u to 0
	for(int i=0;i<N;i++)
		for(int j=i;j<N;j++)
			u[i][j]=A[i][j];//Produce u
	/////////////Produce l///////////
	MAT l(N);							//Delare l
	l=0;								//Set l to 0
	for(int i=0;i<N;i++){
		for(int j=0;j<i;j++)
			l[i][j]=A[i][j];//Produce l
		l[i][i]=1;
	}
	VEC M=fwdSubs(l,b);
	y=bckSubs(u,M);
	return;
}

