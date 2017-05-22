//////VEC.h	 
//////Author:	YUNG-SHAN SU 
//////Date: 3/19 



#ifndef VEC_H
#define VEC_H

class VEC {
	private:
		int dim; // vector length
		double *val; // array to store vector
	public:
		VEC(int n); // uninit constructor, val set to 0
		VEC(const VEC &v1); // copy constructor
		VEC(int n,double *v); // init constructor
		~VEC(); // destructor
		int len(); // dimension of the vector
		VEC &operator-(); // unary operator, negative value
		VEC &operator=(const VEC &v1); // assignment
		VEC &operator=(const double &a); // assignment
		VEC &operator+=(const VEC v1); // V += v1;
		VEC &operator-=(const VEC v1); // V -= v1;
		VEC &operator*=(double a); // V *= dbl;
		VEC &operator/=(double a); // V /= dbl;
		VEC operator+(const VEC &v1); // V + v1
		VEC operator-(const VEC &v1); // V - v1
		double operator*(VEC &v1); // inner product
		VEC operator*(double a); // V * dbl
		VEC operator/(double a); // V / dbl
		VEC absolute();		//absolute(vector)
		double &operator[](int n); // indexing
		void	print();
	friend VEC operator*(double a,const VEC v1); // dbl x V
	friend VEC *newVEC(int n); // create dynamic VEC
	friend VEC absolute(VEC &v1);		//absolute(vector)
};
VEC operator*(double a,const VEC v1);
VEC *newVEC(int n); // create dynamic VEC
//
#endif	
