#ifndef VMATH_HPP
#define VMATH_HPP

#include <cmath>

namespace Vmath {

template<class T>  void Fill( int n, const T alpha,  T *x, const int incx );

template<class T>  void Vmul( int n, const T *x, const int incx, const T *y,
			      const int incy,  T*z, const int incz);

template<class T>  void Smul( int n, const T alpha, const T *x, const int incx,
			      T *y, const int incy);

template<class T>  void Vdiv( int n, const T *x, const int incx, const T *y,
			      const int incy,  T*z, const int incz);

template<class T>  void Sdiv( int n, const T alpha, const T *x,
			      const int incx, T *y, const int incy);

template<class T>  void Vadd( int n, const T *x, const int incx, const T *y,
			      const int incy,  T *z, const int incz);

template<class T>  void Sadd( int n, const T alpha, const T *x,
			      const int incx, T *y, const int incy);

template<class T>  void Vsub( int n, const T *x, const int incx, const T *y,
			      const int incy,  T *z, const int incz);

template<class T>  void Zero(int n, T *x, const int incx);

template<class T>  void Neg( int n, T *x, const int incx);

template<class T> void Vexp(int n, const T *x, const int incx,
			    T *y, const int incy)
{
    while (n--) {
	*y = exp(*x);
	x += incx;
	y += incy;
    }
	
}

template<class T> void Vsqrt(int n, const T *x, const int incx,
			     T *y, const int incy);

template<class T> void Vabs(int n, const T *x, const int incx,
			    T *y, const int incy);

template<class T> void Vvtvp(int n,
			     const T *w, const int incw,
			     const T *x, const int incx,
			     const T *y, const int incy,
			     T *z, const int incz);

template<class T> void Vvtvm(int n, const T *w, const int incw, const T *x,
			     const int incx, const T *y, const int incy,
			     T *z, const int incz);

template<class T> void Svtvp(int n, const T alpha, const T *x,
			     const int incx, const T *y, const int incy,
			     T *z, const int incz);

template<class T> void Svtvm(int n, const T alpha, const T *x,
			     const int incx, const T *y, const int incy,
			     T *z, const int incz);

template<class T> void Vvtvvtp (int n,
				const T* v, int incv,
				const T* w, int incw,
				const T* x, int incx,
				const T* y, int incy,
				T* z, int incz);

template<class T> void Vvtvvtm (int n,
				const T* v, int incv,
				const T* w, int incw,
				const T* x, int incx,
				const T* y, int incy,
				T* z, int incz);

template<class T> void Svtsvtp (int n,
				const T alpha,
				const T* x, int incx,
				const T beta,
				const T* y, int incy,
				T* z, int incz);

template<class T> void Vstvpp(int n,
			      const T alpha,
			      const T* v, int incv,
			      const T* w, int incw,
			      const T* x, int incx,
			      T* z, int incz);

template<typename T> void Vcopy(int n, const T*x, const int incx,
				T *y, const int incy);

template<class T> void Vlog(int n, const T *x, const int incx,
				T *y, const int incy)
{
    while (n--)
	{
	    *y = log( *x );
	    x += incx;
	    y += incy;
	}
}

template<class T> void Vpow(int n, const T *x, const int incx,
                const T f, T *y, const int incy)
{
    while (n--)
        {
            *y = pow( *x, f );
            x += incx;
            y += incy;
        }
}    

}
#endif //VMATH_HPP
