#include <cstring>
#include <cmath>
#include "SharedArray.hpp"
#include "Vmath.hpp"

namespace Vmath {

template<class T>  void Fill( int n, const T alpha,  T *x, const int incx )
{
    while( n-- ) {
	*x = alpha;
	x += incx;
    }
}

template void Fill( int n, const double alpha,  double *x, const int incx );
template void Fill( int n, const float alpha,  float *x, const int incx );

template<class T>  void Vmul( int n, const T *x, const int incx, const T *y,
			      const int incy,  T*z, const int incz)
{
    ++n;
    if (incx == 1 && incy == 1 && incz == 1) {
	while( --n )
	    {
		*z = (*x) * (*y);
		++x;
		++y;
		++z;
	    }
    }
    else {
	while( --n )
	    {
		*z = (*x) * (*y);
		x += incx;
		y += incy;
		z += incz;
	    }
    }
}

template void Vmul( int n, const double *x, const int incx, const double *y,
			      const int incy,  double*z, const int incz);
template void Vmul( int n, const float *x, const int incx, const float *y,
			      const int incy,  float*z, const int incz);

template<class T>  void Smul( int n, const T alpha, const T *x, const int incx,
			      T *y, const int incy)
{
    ++n;
    if (incx == 1 && incy == 1) {
	while( --n )
	    {
		*y = alpha * (*x);
		++x;
		++y;
	    }
    }
    else {

	while( --n )
	    {
		*y = alpha * (*x);
		x += incx;
		y += incy;
	    }
    }
}

template void Smul( int n, const double alpha, const double *x, const int incx,
		    double *y, const int incy);
template void Smul( int n, const float alpha, const float *x, const int incx,
		    float *y, const int incy);

template<class T>  void Vdiv( int n, const T *x, const int incx, const T *y,
			      const int incy,  T*z, const int incz)
{
    ++n;
    if (incx == 1 && incy == 1) {
	while( --n )
	    {
		*z = (*x) / (*y);
		++x;
		++y;
		++z;
	    }
    }
    else {
	while( --n )
	    {
		*z = (*x) / (*y);
		x += incx;
		y += incy;
		z += incz;
	    }
    }
}

template void Vdiv( int n, const double *x, const int incx, const double *y,
		    const int incy,  double*z, const int incz);
template void Vdiv( int n, const float *x, const int incx, const float *y,
		    const int incy,  float*z, const int incz);


template<class T>  void Sdiv( int n, const T alpha, const T *x,
			      const int incx, T *y, const int incy)
{
    ++n;
    if (incx == 1 && incy == 1) {
	while( --n ) {
	    *y = alpha / (*x);
	    ++x;
	    ++y;
	}
    }
    else {
	while( --n ) {
	    *y = alpha / (*x);
	    x += incx;
	    y += incy;
	}
    }
}

template void Sdiv( int n, const double alpha, const double *x,
		    const int incx, double *y, const int incy);
template void Sdiv( int n, const float alpha, const float *x,
		    const int incx, float *y, const int incy);

template<class T>  void Vadd( int n, const T *x, const int incx, const T *y,
			      const int incy,  T *z, const int incz)
{
    while( n-- ) {
	*z = (*x) + (*y);
	x += incx;
	y += incy;
	z += incz;
    }
}

template void Vadd( int n, const double *x, const int incx, const double *y,
		    const int incy,  double *z, const int incz);
template void Vadd( int n, const float *x, const int incx, const float *y,
		    const int incy,  float *z, const int incz);

template<class T>  void Sadd( int n, const T alpha, const T *x,
			      const int incx, T *y, const int incy)
{
    ++n;
    if (incx == 1 && incy == 1) {
	while( --n ) {
	    *y = alpha + (*x);
	    ++x;
	    ++y;
	}
    }
    else {

	while( --n ) {
	    *y = alpha + (*x);
	    x += incx;
	    y += incy;
	}
    }
}

template void Sadd( int n, const double alpha, const double *x,
		    const int incx, double *y, const int incy);
template void Sadd( int n, const float alpha, const float *x,
		    const int incx, float *y, const int incy);

template<class T>  void Vsub( int n, const T *x, const int incx, const T *y,
			      const int incy,  T *z, const int incz)
{
    ++n;
    if (incx == 1 && incy == 1 && incz == 1) {
	while( --n ) {
	    *z = (*x) - (*y);
	    ++x;
	    ++y;
	    ++z;
	}
    }
    else {
	while( --n ) {
	    *z = (*x) - (*y);
	    x += incx;
	    y += incy;
	    z += incz;
	}
    }
}

template void Vsub( int n, const double *x, const int incx, const double *y,
			      const int incy,  double *z, const int incz);
template void Vsub( int n, const float *x, const int incx, const float *y,
			      const int incy,  float *z, const int incz);

template<class T>  void Zero(int n, T *x, const int incx)
{
    if(incx == 1) {
	
	std::memset(x,'\0', n*sizeof(T));
    }
    else {
	T zero = 0;
	++n;
	while(--n)
            {
		*x = zero;
		x+=incx;
            }
    }
}

template void Zero(int n, double *x, const int incx);
template void Zero(int n, float *x, const int incx); 

template<class T>  void Neg( int n, T *x, const int incx)
{
    while( n-- ) {
	*x = -(*x);
	x += incx;
    }
}

template void Neg( int n, double *x, const int incx);
template void Neg( int n, float *x, const int incx);

    /// \brief sqrt y = sqrt(x)
template<class T> void Vsqrt(int n, const T *x, const int incx,
			     T *y, const int incy)
{
    while (n--) {
	*y  = sqrt( *x );
	x  += incx;
	y  += incy;
    }
}

template void Vsqrt(int n, const double *x, const int incx,
			     double *y, const int incy);
template void Vsqrt(int n, const float *x, const int incx,
			     float *y, const int incy);

    /// \brief vabs: y = |x|
template<class T> void Vabs(int n, const T *x, const int incx,
			    T *y, const int incy)
{
    while( n-- ) {
	*y = ( *x >0)? *x:-(*x);
	x += incx;
	y += incy;
    }
}

template void Vabs(int n, const double *x, const int incx,
		   double *y, const int incy);
template void Vabs(int n, const float *x, const int incx,
		   float *y, const int incy);

    /********** Triad  routines  ***********************/

    /// \brief  vvtvp (vector times vector plus vector): z = w*x + y
template<class T> void Vvtvp(int n,
			     const T *w, const int incw,
			     const T *x, const int incx,
			     const T *y, const int incy,
			     T *z, const int incz)
{
    while( n-- ) {
	*z = (*w) * (*x) + (*y);
	w += incw;
	x += incx;
	y += incy;
	z += incz;
    }
}


template void Vvtvp(int n,
		    const double *w, const int incw,
		    const double *x, const int incx,
		    const double *y, const int incy,
		    double *z, const int incz);
template void Vvtvp(int n,
		    const float *w, const int incw,
		    const float *x, const int incx,
		    const float *y, const int incy,
		    float *z, const int incz);

    /// \brief vvtvm (vector times vector plus vector): z = w*x - y
template<class T> void Vvtvm(int n, const T *w, const int incw, const T *x,
			     const int incx, const T *y, const int incy,
			     T *z, const int incz)
{
    while( n-- ) {
	*z = (*w) * (*x) - (*y);
	w += incw;
	x += incx;
	y += incy;
	z += incz;
    }
}

template void Vvtvm(int n, const double *w, const int incw, const double *x,
			     const int incx, const double *y, const int incy,
		    double *z, const int incz);
template void Vvtvm(int n, const float *w, const int incw, const float *x,
			     const int incx, const float *y, const int incy,
		    float *z, const int incz);

    /// \brief  svtvp (scalar times vector plus vector): z = alpha*x + y
template<class T> void Svtvp(int n, const T alpha, const T *x,
			     const int incx, const T *y, const int incy,
			     T *z, const int incz)
{
    ++n;
    if (incx == 1 && incy == 1 && incz == 1) {
	while( --n ) {
	    *z = alpha * (*x) + (*y);
	    ++x;
	    ++y;
	    ++z;
	}
    }
    else {
	while( --n ) {
	    *z = alpha * (*x) + (*y);
	    x += incx;
	    y += incy;
	    z += incz;
	}
    }
}

template void Svtvp(int n, const double alpha, const double *x,
		    const int incx, const double *y, const int incy,
		    double *z, const int incz);
template void Svtvp(int n, const float alpha, const float *x,
		    const int incx, const float *y, const int incy,
		    float *z, const int incz);

    /// \brief  svtvp (scalar times vector plus vector): z = alpha*x - y
template<class T> void Svtvm(int n, const T alpha, const T *x,
			     const int incx, const T *y, const int incy,
			     T *z, const int incz)
{
    while( n-- ) {
	*z = alpha * (*x) - (*y);
	x += incx;
	y += incy;
	z += incz;
    }
}

template void Svtvm(int n, const double alpha, const double *x,
		    const int incx, const double *y, const int incy,
		    double *z, const int incz);
template void Svtvm(int n, const float alpha, const float *x,
		    const int incx, const float *y, const int incy,
		    float *z, const int incz);

    /// \brief  vvtvvtp (vector times vector plus vector times vector):
    // z = v*w + x*y
template<class T> void Vvtvvtp (int n,
				const T* v, int incv,
				const T* w, int incw,
				const T* x, int incx,
				const T* y, int incy,
				T* z, int incz)
{
    while( n-- ){
	*z = (*v) * (*w) + (*x) * (*y);
	v += incv;
	w += incw;
	x += incx;
	y += incy;
	z += incz;
    }
}

template void Vvtvvtp (int n,
		       const double* v, int incv,
		       const double* w, int incw,
		       const double* x, int incx,
		       const double* y, int incy,
		       double* z, int incz);
template void Vvtvvtp (int n,
		       const float* v, int incv,
		       const float* w, int incw,
		       const float* x, int incx,
		       const float* y, int incy,
		       float* z, int incz);

    /// \brief  vvtvvtm (vector times vector minus vector times vector):
    // z = v*w - x*y
template<class T> void Vvtvvtm (int n,
				const T* v, int incv,
				const T* w, int incw,
				const T* x, int incx,
				const T* y, int incy,
				T* z, int incz)
{
    while( n-- ) {

	*z = (*v) * (*w) - (*x) * (*y);
	v += incv;
	w += incw;
	x += incx;
	y += incy;
	z += incz;
    }
}

template void Vvtvvtm (int n,
		       const double* v, int incv,
		       const double* w, int incw,
		       const double* x, int incx,
		       const double* y, int incy,
		       double* z, int incz);
template void Vvtvvtm (int n,
		       const float* v, int incv,
		       const float* w, int incw,
		       const float* x, int incx,
		       const float* y, int incy,
		       float* z, int incz);

    /// \brief  vvtvvtp (scalar times vector plus scalar times vector):
    // z = alpha*x + beta*y
template<class T> void Svtsvtp (int n,
				const T alpha,
				const T* x, int incx,
				const T beta,
				const T* y, int incy,
				T* z, int incz)
{
    while( n-- ) {

	*z = alpha * (*x) + beta * (*y);
	x += incx;
	y += incy;
	z += incz;
    }
}

template void Svtsvtp (int n,
		       const double alpha,
		       const double* x, int incx,
		       const double beta,
		       const double* y, int incy,
		       double* z, int incz);
template void Svtsvtp (int n,
		       const float alpha,
		       const float* x, int incx,
		       const float beta,
		       const float* y, int incy,
		       float* z, int incz);

    /// \brief  Vstvpp (scalar times vector plus vector plus vector):
    // z = v*w + x*y
template<class T> void Vstvpp(int n,
			      const T alpha,
			      const T* v, int incv,
			      const T* w, int incw,
			      const T* x, int incx,
			      T* z, int incz)
{
    while( n-- ) {
	*z = alpha * (*v) + (*w) + (*x);
	v += incv;
	w += incw;
	x += incx;
	z += incz;
    }
}

template void Vstvpp(int n, const double alpha, const double* v, int
		     incv, const double* w, int incw, const double* x,
		     int incx, double* z, int incz);
template void Vstvpp(int n, const float alpha, const float* v, int
		     incv, const float* w, int incw, const float* x,
		     int incx, float* z, int incz);

template<typename T> void Vcopy(int n, const T *x, const int incx, T
*y, const int incy) {
    if( incx ==1 && incy == 1) {
	memcpy(y,x,n*sizeof(T));
    }
    else {
	while( n-- ) {
	    *y = *x;
	    x += incx;
	    y += incy;
	}
    }
}

template void Vcopy(int n, const double *x, const int incx, double *y,
const int incy);
template void Vcopy(int n, const float *x, const int incx, float *y,
const int incy);
}
