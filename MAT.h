//////MAT.cpp		 
//////Author:	YUNG-SHAN SU 
//////Date: 3/25 




// matrix class
#ifndef MAT_H
#define MAT_H
#include "VEC.h"
class MAT {
	private:
		int n; // define nxn matrix
		int ROW;
		int COLUMN;
		VEC **va; // array of n pointers to vectors
	public:
		MAT(int dim); // uninit constructor
		MAT(int row,int column); // uninit constructor
		MAT(const MAT &m1); // copy constructor
		MAT(int dim,double *v); // init constructor
		~MAT(); // destructor
		int dim(); // return dimension of the matrix
		int getrow();
        int getcolumn();
		MAT tpose(); // transpose
		MAT &operator-(); // unary operator, negative value
		MAT &operator=(MAT m1); // assignment
		MAT &operator=(double a); // assignment
		MAT &operator+=(MAT &m1); // m += m1;
		MAT &operator-=(MAT &m1); // m -= m1;
		MAT &operator*=(double a); // m *= dbl;
		MAT &operator/=(double a); // m /= dbl;
		MAT operator+(MAT m1); // m1 + m2
		MAT operator-(MAT m1); // m1 - m2
		MAT operator*(MAT m1); // m1 * m2
		VEC & operator[](int m); // m'th row
		VEC operator*(VEC &v1); // m x v1
		MAT operator*(double a); // m * dbl
		MAT operator/(double a); // m / dbl
		void print();
	friend MAT operator*(double a,MAT &m1); // dbl x m
	friend VEC operator*(VEC &v1,MAT &m1); // vT x m
};
MAT operator*(double a,const MAT &m1); // dbl x m
VEC operator*(VEC &v1,MAT &m1); // vT x m


/////////HW2 function////////////
MAT &luFact(MAT &m1); // LU decomposition
VEC fwdSubs(MAT &m1,VEC b); // forward substitution
VEC bckSubs(MAT &m1,VEC b); // backward substitution
/////////////////////////////////

/////////HW3 function////////////
MAT  getResistorNetworksMatrix(int dim); //Get the linear system network
/////////////////////////////////

////////HW4 function////////////
int jacobi(MAT &A,VEC b,VEC &x,int maxIter,double tol);
int gaussSeidel(MAT &A,VEC b,VEC &x,int maxIter,double tol);
int sgs(MAT &A,VEC b,VEC &x,int maxIter,double tol);
double one_norm(VEC &va);
double two_norm(VEC &va);
double max_norm(VEC &va);

////////HW5 function////////////
int cg(MAT &A,VEC b,VEC &x,int maxIter, double tol);


/////////HW6 function//////////////
int EVpwr(MAT &A,VEC &q0,double &lambda,double tol,int maxiter);
int EViPwr(MAT &A,VEC &q0,double &lambda,double tol,int maxiter);
int EViPwrShft(MAT &A,VEC &q0,double &lambda,double mu,double tol,int maxiter);
int QR_shift(MAT &A,VEC &lambda,double &condition_number,double mu,double tol,int maxiter);

double epson_1(double &lambda1,double &lambda2);
double epson_2(VEC &q1,VEC &q2);
double epson_3(MAT &A,VEC &q,double &lambda);
double epson_4(MAT &A,VEC &q,double &lambda);

///////HW7 function//////////////
int EVqr(MAT &A,double tol,int maxiter);
int EVqrShifted(MAT &A,double mu,double tol,int maxiter);

#endif

///////HW8 function//////////////
double NEV(double x,VEC &XS,VEC &YS,int n);		//Non-recursive Neville’s Algorithm
double Lagrange(double x,VEC &XDATA,VEC &YDATA);//Lagrange Algorithm
////////HW9 function/////////////////
void splineM(int N,VEC &X,VEC &Y,VEC &M); // generate spline momentum M
double spline(double x,int N,VEC &X,VEC &Y,VEC &M); // spline interp at x



////////HW11 finction//////////////////
void getLinearSolution(MAT A,VEC & y,VEC b); 

////////HW12 function
VEC function12_forwardEuler(VEC &x_t, double step);
VEC function12_backwardEuler(VEC &x_t,double step);
VEC function12_trep(VEC &x_t,double step);
MAT forwardEuler(VEC &initial,double start, double end, double step);
MAT backwardEuler(VEC &initial,double start, double end, double step);
MAT trep(VEC &initial,double start, double end, double step);
VEC columnMax(MAT &m1);
VEC columnMin(MAT &m1);