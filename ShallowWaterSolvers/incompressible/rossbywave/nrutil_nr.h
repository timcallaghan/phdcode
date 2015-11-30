#ifndef _NR_UTIL_H_
#define _NR_UTIL_H_

#include <string>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
using namespace std;

typedef double DP;

template<class T>
inline const T SQR(const T a) {return a*a;}

template<class T>
inline const T MAX(const T &a, const T &b)
        {return b > a ? (b) : (a);}

inline float MAX(const double &a, const float &b)
        {return b > a ? (b) : float(a);}

inline float MAX(const float &a, const double &b)
        {return b > a ? float(b) : (a);}

template<class T>
inline const T MIN(const T &a, const T &b)
        {return b < a ? (b) : (a);}

inline float MIN(const double &a, const float &b)
        {return b < a ? (b) : float(a);}

inline float MIN(const float &a, const double &b)
        {return b < a ? float(b) : (a);}

template<class T>
inline const T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

template<class T>
inline void SWAP(T &a, T &b)
	{T dum=a; a=b; b=dum;}

namespace NR {
	inline void nrerror(const string error_text)
	// Numerical Recipes standard error handler
	{
		cerr << "Numerical Recipes run-time error..." << endl;
		cerr << error_text << endl;
		cerr << "...now exiting to system..." << endl;
		exit(1);
	}
}

namespace NR {
	inline void nrwarn(const string error_text)
	// Numerical Recipes standard error handler
	{
		cerr << "Numerical Recipes run-time warning..." << endl;
		cerr << error_text << endl;
	}
}

template <class T>
class NRVec {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	NRVec();
	explicit NRVec(int n);		// Zero-based array
	NRVec(const T &a, int n);	//initialize to constant value
	NRVec(const T *a, int n);	// Initialize to array
	NRVec(const NRVec &rhs);	// Copy constructor
	NRVec & operator=(const NRVec &rhs);	//assignment
	NRVec & operator=(const T &a);	//assign a to every element
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	~NRVec();
};

template <class T>
NRVec<T>::NRVec() : nn(0), v(0) {}

template <class T>
NRVec<T>::NRVec(int n) : nn(n), v(new T[n]) {}

template <class T>
NRVec<T>::NRVec(const T& a, int n) : nn(n), v(new T[n])
{
	for(int i=0; i<n; i++)
		v[i] = a;
}

template <class T>
NRVec<T>::NRVec(const T *a, int n) : nn(n), v(new T[n])
{
	for(int i=0; i<n; i++)
		v[i] = *a++;
}

template <class T>
NRVec<T>::NRVec(const NRVec<T> &rhs) : nn(rhs.nn), v(new T[nn])
{
	for(int i=0; i<nn; i++)
		v[i] = rhs[i];
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const NRVec<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.nn) {
			if (v != 0) delete [] (v);
			nn=rhs.nn;
			v= new T[nn];
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i<nn; i++)
		v[i]=a;
	return *this;
}

template <class T>
inline T & NRVec<T>::operator[](const int i)	//subscripting
{
	return v[i];
}

template <class T>
inline const T & NRVec<T>::operator[](const int i) const	//subscripting
{
	return v[i];
}

template <class T>
inline int NRVec<T>::size() const
{
	return nn;
}

template <class T>
NRVec<T>::~NRVec()
{
	if (v != 0)
		delete[] (v);
}

/* My own overloaded opperators and functions for working
with vectors*/

template <class T>
NRVec<T> operator+(const NRVec<T> &A, const NRVec<T> &B)
// Overloaded addition
{
    if (A.size() != B.size())
	{
		NR::nrerror("Can't add two vectors of unequal length");
	}

	int N=A.size();
    NRVec<T> tmp(N);
    for (int i=0; i<N; i++)
            tmp[i] = A[i] + B[i];

    return tmp;
}

template <class T>
NRVec<T> operator-(const NRVec<T> &A, const NRVec<T> &B)
// Overloaded subtraction
{
    if (A.size() != B.size())
	{
		NR::nrerror("Can't subtract two vectors of unequal length");
	}

	int N=A.size();
    NRVec<T> tmp(N);
    for (int i=0; i<N; i++)
            tmp[i] = A[i] - B[i];

    return tmp;
}

template <class T>
T operator*(const NRVec<T> &A, const NRVec<T> &B)
// Overloaded multiplication for vector dot product
{
	if (A.size() != B.size())
	{
		NR::nrerror("Can't dot product two vectors of unequal length");
	}
    
	int N = A.size();
	T sum=0.0;
	for (int i=0;i<N;i++)
	{
		sum += A[i]*B[i];
	}
	return sum;
}

template <class T>
NRVec<T> operator*(const NRVec<T> &A, const T &a) 
// Post multiply by a constant.
{
	int N = A.size();
	NRVec<T> tmp(N);
	for (int i=0;i<N;i++)
	{
		tmp[i]=A[i]*a;
	}
	return tmp;

}

template <class T>
NRVec<T> operator*(const NRVec<T> &A, const int &a) 
// Post multiply by an integer.
{
	int N = A.size();
	NRVec<T> tmp(N);
	for (int i=0;i<N;i++)
	{
		tmp[i]=A[i]*a;
	}
	return tmp;

}

template <class T>
NRVec<T> operator*(const T &a, const NRVec<T> &A) 
// Pre multiply by a constant.
{
	int N = A.size();
	NRVec<T> tmp(N);
	for (int i=0;i<N;i++)
	{
		tmp[i]=A[i]*a;
	}
	return tmp;

}

template <class T>
NRVec<T> operator*(const int &a, const NRVec<T> &A) 
// Pre multiply by an integer.
{
	int N = A.size();
	NRVec<T> tmp(N);
	for (int i=0;i<N;i++)
	{
		tmp[i]=A[i]*a;
	}
	return tmp;

}


/* ***************************  I/O  ********************************/

template <class T>
std::ostream& operator<<(std::ostream &s, const NRVec<T> &A)
{
    int N=A.size();

    for (int i=0; i<N; i++)
        s   << A[i] << " " << endl;
    s << endl;

    return s;
}

// End of my own overloaded opperators

template <class T>
class NRMat {
private:
	int nn;
	int mm;
	T **v;
public:
	NRMat();
	NRMat(int n, int m);			// Zero-based array
	NRMat(const T &a, int n, int m);	//Initialize to constant
	NRMat(const T *a, int n, int m);	// Initialize to array
	NRMat(const NRMat &rhs);		// Copy constructor
	NRMat & operator=(const NRMat &rhs);	//assignment
	NRMat & operator=(const T &a);		//assign a to every element
	inline T* operator[](const int i);	//subscripting: pointer to row i
	inline const T* operator[](const int i) const;
	inline int nrows() const;
	inline int ncols() const;
	~NRMat();
};

template <class T>
NRMat<T>::NRMat() : nn(0), mm(0), v(0) {}

template <class T>
NRMat<T>::NRMat(int n, int m) : nn(n), mm(m), v(new T*[n])
{
	v[0] = new T[m*n];
	for (int i=1; i< n; i++)
		v[i] = v[i-1] + m;
}

template <class T>
NRMat<T>::NRMat(const T &a, int n, int m) : nn(n), mm(m), v(new T*[n])
{
	int i,j;
	v[0] = new T[m*n];
	for (i=1; i< n; i++)
		v[i] = v[i-1] + m;
	for (i=0; i< n; i++)
		for (j=0; j<m; j++)
			v[i][j] = a;
}

template <class T>
NRMat<T>::NRMat(const T *a, int n, int m) : nn(n), mm(m), v(new T*[n])
{
	int i,j;
	v[0] = new T[m*n];
	for (i=1; i< n; i++)
		v[i] = v[i-1] + m;
	for (i=0; i< n; i++)
		for (j=0; j<m; j++)
			v[i][j] = *a++;
}

template <class T>
NRMat<T>::NRMat(const NRMat &rhs) : nn(rhs.nn), mm(rhs.mm), v(new T*[nn])
{
	int i,j;
	v[0] = new T[mm*nn];
	for (i=1; i< nn; i++)
		v[i] = v[i-1] + mm;
	for (i=0; i< nn; i++)
		for (j=0; j<mm; j++)
			v[i][j] = rhs[i][j];
}

template <class T>
NRMat<T> & NRMat<T>::operator=(const NRMat<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		int i,j;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (v != 0) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nn;
			mm=rhs.mm;
			v = new T*[nn];
			v[0] = new T[mm*nn];
		}
		for (i=1; i< nn; i++)
			v[i] = v[i-1] + mm;
		for (i=0; i< nn; i++)
			for (j=0; j<mm; j++)
				v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
NRMat<T> & NRMat<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i< nn; i++)
		for (int j=0; j<mm; j++)
			v[i][j] = a;
	return *this;
}

template <class T>
inline T* NRMat<T>::operator[](const int i)	//subscripting: pointer to row i
{
	return v[i];
}

template <class T>
inline const T* NRMat<T>::operator[](const int i) const
{
	return v[i];
}

template <class T>
inline int NRMat<T>::nrows() const
{
	return nn;
}

template <class T>
inline int NRMat<T>::ncols() const
{
	return mm;
}

template <class T>
NRMat<T>::~NRMat()
{
	if (v != 0) {
		delete[] (v[0]);
		delete[] (v);
	}
}

// My own overloaded operators and functions for
// easy matrix and vector manipulation.

template <class T>
NRMat<T> operator+(const NRMat<T> &A, const NRMat<T> &B)
// overloaded matrix addition
{
    int M = A.nrows();
    int N = A.ncols();

	if ( (M != B.nrows()) || (N != B.ncols()) )
	{
		NR::nrerror("Can't add two matricies of unequal size");
	}
    
    NRMat<T> tmp(M,N);
    int i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] = A[i][j] + B[i][j];

    return tmp;
}

template <class T>
NRMat<T> operator-(const NRMat<T> &A, const NRMat<T> &B)
// overloaded matrix subtraction
{
    int M = A.nrows();
    int N = A.ncols();

	if ( (M != B.nrows()) || (N != B.ncols()) )
	{
		NR::nrerror("Can't subtract two matricies of unequal size");
	}
    
    NRMat<T> tmp(M,N);
    int i,j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            tmp[i][j] = A[i][j] - B[i][j];

    return tmp;
}

template <class T>
NRMat<T> transposeC(const NRMat<T> &A)
// Matrix transpose...create a new one without altering A
{
    int M = A.nrows();
    int N = A.ncols();

    NRMat<T> S(N,M);
    int i, j;

    for (i=0; i<M; i++)
        for (j=0; j<N; j++)
            S[j][i] = A[i][j];

    return S;
}

template <class T>
void transposeR(NRMat<T> &A)
// Matrix transpose...replace A with it's transpoe
{
    int M = A.nrows();
    int N = A.ncols();
	// Indexes used in the tranpose operation
	int i,j;
	// Test if we can do quick transpose...only possible
	// if M=N. If it's not, we have to return an error.
	if (M != N)
	{
		NR::nrerror("Can't take replacement tranpose...will require resize...");
	}
	// Temp value used to hold the reflection of an element
	T temp;
	for (i=0;i<M;i++)
	{
		for (j=i;j<N;j++)
		{
			// Store current element
			temp=A[i][j];
			// Reflect element and store
			A[i][j]=A[j][i];
			// Store the temp in the reflection
			A[j][i]=temp;
		}
	}
}

template <class T>
inline NRMat<T> matmult(const NRMat<T>  &A, const NRMat<T> &B)
// Matrix multiplication
{

	if (A.ncols() != B.nrows())
	{
		NR::nrerror("Can't multiply matrices with different numbers of columns and rows");
	}

    int M = A.nrows();
    int N = A.ncols();
    int K = B.ncols();

    NRMat<T> tmp(M,K);
    T sum;

    for (int i=0; i<M; i++)
    for (int k=0; k<K; k++)
    {
        sum = 0;
        for (int j=0; j<N; j++)
            sum = sum +  A[i][j] * B[j][k];

        tmp[i][k] = sum; 
    }

    return tmp;
}

template <class T>
inline NRMat<T> operator*(const NRMat<T>  &A, const NRMat<T> &B)
// Overloaded matrix multiplication
{
    return matmult(A,B);
}

template <class T>
NRVec<T> matmult(const NRMat<T>  &A, const NRVec<T> &x)
// Matrix vector multiplication
{

	if (A.ncols() != x.size())
	{
		NR::nrerror("Matrix vector multiplication dimensions must agree");
	}

    int M = A.nrows();
    int N = A.ncols();

    NRVec<T> tmp(M);
    T sum;

    for (int i=0; i<M; i++)
    {
        sum = 0;
        const T* rowi = A[i];
        for (int j=0; j<N; j++)
            sum = sum +  rowi[j] * x[j];

        tmp[i] = sum; 
    }

    return tmp;
}

template <class T>
inline NRVec<T> operator*(const NRMat<T>  &A, const NRVec<T> &x)
// Overloaded matrix vector multiplication
{
    return matmult(A,x);
}

/* ***************************  I/O  ********************************/

template <class T>
std::ostream& operator<<(std::ostream &s, const NRMat<T> &A)
{
    int M=A.nrows();
    int N=A.ncols();

    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
        {
            s << setw(12) << A[i][j] << " ";
        }
        s << endl;
    }


    return s;
}

template <class T>
class NRMat3d {
private:
	int nn;
	int mm;
	int kk;
	T ***v;
public:
	NRMat3d();
	NRMat3d(int n, int m, int k);
	inline T** operator[](const int i);	//subscripting: pointer to row i
	inline const T* const * operator[](const int i) const;
	inline int dim1() const;
	inline int dim2() const;
	inline int dim3() const;
	~NRMat3d();
};

template <class T>
NRMat3d<T>::NRMat3d(): nn(0), mm(0), kk(0), v(0) {}

template <class T>
NRMat3d<T>::NRMat3d(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
	int i,j;
	v[0] = new T*[n*m];
	v[0][0] = new T[n*m*k];
	for(j=1; j<m; j++)
		v[0][j] = v[0][j-1] + k;
	for(i=1; i<n; i++) {
		v[i] = v[i-1] + m;
		v[i][0] = v[i-1][0] + m*k;
		for(j=1; j<m; j++)
			v[i][j] = v[i][j-1] + k;
	}
}

template <class T>
inline T** NRMat3d<T>::operator[](const int i) //subscripting: pointer to row i
{
	return v[i];
}

template <class T>
inline const T* const * NRMat3d<T>::operator[](const int i) const
{
	return v[i];
}

template <class T>
inline int NRMat3d<T>::dim1() const
{
	return nn;
}

template <class T>
inline int NRMat3d<T>::dim2() const
{
	return mm;
}

template <class T>
inline int NRMat3d<T>::dim3() const
{
	return kk;
}

template <class T>
NRMat3d<T>::~NRMat3d()
{
	if (v != 0) {
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	}
}

/* ***************************  I/O  ********************************/
// template for easy stream output of the 3d Matrix
template <class T>
std::ostream& operator<<(std::ostream &s, const NRMat3d<T> &A)
{
    int M=A.dim1();
    int N=A.dim2();
	int O=A.dim3();

	for (int k=0;k<O;k++)
	{
		for (int i=0; i<M; i++)
		{
			for (int j=0; j<N; j++)
			{
				s << setw(12) << A[i][j][k] << " ";
			}
			s << endl;
		}
		s << endl;
	}
    return s;
}

//The next 3 classes are used in artihmetic coding, Huffman coding, and
//wavelet transforms respectively. This is as good a place as any to put them!

class arithcode {
private:
	NRVec<unsigned long> *ilob_p,*iupb_p,*ncumfq_p;
public:
	NRVec<unsigned long> &ilob,&iupb,&ncumfq;
	unsigned long jdif,nc,minint,nch,ncum,nrad;
	arithcode(unsigned long n1, unsigned long n2, unsigned long n3)
		: ilob_p(new NRVec<unsigned long>(n1)),
		iupb_p(new NRVec<unsigned long>(n2)),
		ncumfq_p(new NRVec<unsigned long>(n3)),
		ilob(*ilob_p),iupb(*iupb_p),ncumfq(*ncumfq_p) {}
	~arithcode() {
		if (ilob_p != 0) delete ilob_p;
		if (iupb_p != 0) delete iupb_p;
		if (ncumfq_p != 0) delete ncumfq_p;
	}
};

class huffcode {
private:
	NRVec<unsigned long> *icod_p,*ncod_p,*left_p,*right_p;
public:
	NRVec<unsigned long> &icod,&ncod,&left,&right;
	int nch,nodemax;
	huffcode(unsigned long n1, unsigned long n2, unsigned long n3,
		unsigned long n4) :
		icod_p(new NRVec<unsigned long>(n1)),
		ncod_p(new NRVec<unsigned long>(n2)),
		left_p(new NRVec<unsigned long>(n3)),
		right_p(new NRVec<unsigned long>(n4)),
		icod(*icod_p),ncod(*ncod_p),left(*left_p),right(*right_p) {}
	~huffcode() {
		if (icod_p != 0) delete icod_p;
		if (ncod_p != 0) delete ncod_p;
		if (left_p != 0) delete left_p;
		if (right_p != 0) delete right_p;
	}
};

class wavefilt {
private:
	NRVec<DP> *cc_p,*cr_p;
public:
	int ncof,ioff,joff;
	NRVec<DP> &cc,&cr;
	wavefilt() : cc(*cc_p),cr(*cr_p) {}
	wavefilt(const DP *a, const int n) :  //initialize to array
		cc_p(new NRVec<DP>(n)),cr_p(new NRVec<DP>(n)),
		ncof(n),ioff(-(n >> 1)),joff(-(n >> 1)),cc(*cc_p),cr(*cr_p) {
			int i;
			for (i=0; i<n; i++)
				cc[i] = *a++;
			DP sig = -1.0;
			for (i=0; i<n; i++) {
				cr[n-1-i]=sig*cc[i];
				sig = -sig;
			}
	}
	~wavefilt() {
		if (cc_p != 0) delete cc_p;
		if (cr_p != 0) delete cr_p;
	}
};

//Overloaded complex operations to handle mixed float and double
//This takes care of e.g. 1.0/z, z complex<float>

inline const complex<float> operator+(const double &a,
	const complex<float> &b) { return float(a)+b; }

inline const complex<float> operator+(const complex<float> &a,
	const double &b) { return a+float(b); }

inline const complex<float> operator-(const double &a,
	const complex<float> &b) { return float(a)-b; }

inline const complex<float> operator-(const complex<float> &a,
	const double &b) { return a-float(b); }

inline const complex<float> operator*(const double &a,
	const complex<float> &b) { return float(a)*b; }

inline const complex<float> operator*(const complex<float> &a,
	const double &b) { return a*float(b); }

inline const complex<float> operator/(const double &a,
	const complex<float> &b) { return float(a)/b; }

inline const complex<float> operator/(const complex<float> &a,
	const double &b) { return a/float(b); }

//some compilers choke on pow(float,double) in single precision. also atan2

inline float pow (float x, double y) {return pow(double(x),y);}
inline float pow (double x, float y) {return pow(x,double(y));}
inline float atan2 (float x, double y) {return atan2(double(x),y);}
inline float atan2 (double x, float y) {return atan2(x,double(y));}
#endif /* _NR_UTIL_H_ */


