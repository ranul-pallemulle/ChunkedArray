#ifndef VMATHARRAY_HPP
#define VMATHARRAY_HPP

#include "Vmath.hpp"
#include "SharedArray.hpp"

namespace Vmath
{

template<class T>  void Fill( int n, const T alpha,  Array<OneD, T> &x, const int incx ) {
    Fill(n,alpha,&x[0],incx);
}

template<class T>  void Vmul( int n, const Array<OneD, const T> &x, const int incx, const Array<OneD, const T> &y, const int incy,  Array<OneD,T> &z, const int incz) {
    Vmul(n,&x[0],incx,&y[0],incy,&z[0],incz);
}

template<class T>  void Smul( int n, const T alpha, const Array<OneD,const T> &x,  const int incx,  Array<OneD,T>  &y, const int incy) {
    Smul(n,alpha, &x[0],incx,&y[0],incy);
}

template<class T>  void Vdiv( int n, const Array<OneD,const T> &x, const int incx, const Array<OneD,const T> &y, const int incy,  Array<OneD,T> &z, const int incz) {
    Vdiv(n,&x[0],incx,&y[0],incy,&z[0],incz);
}
        
        /// \brief Scalar multiply  y = alpha/y
template<class T>  void Sdiv( int n, const T alpha, const Array<OneD,const T> &x, const int incx,  Array<OneD,T> &y, const int incy) {
    Sdiv(n,alpha,&x[0],incx,&y[0],incy);
}
        
        /// \brief Add vector z = x+y
template<class T>  void Vadd( int n, const Array<OneD,const T> &x, const int incx, const Array<OneD,const T> &y,  const int incy,  Array<OneD,T> &z, const int incz) {
    Vadd(n,&x[0],incx,&y[0],incy,&z[0],incz);
}
    
        /// \brief Add vector y = alpha + x
template<class T>  void Sadd( int n, const T alpha, const Array<OneD,const T> &x,const int incx, Array<OneD,T> &y, const int incy)
{
    Sadd(n,alpha,&x[0],incx,&y[0],incy);
}
    
        /// \brief Subtract vector z = x-y
template<class T>  void Vsub( int n, const Array<OneD,const T> &x, const int incx, const Array<OneD,const T> &y, const int incy,  Array<OneD,T> &z, const int incz)
{
    Vsub(n,&x[0],incx,&y[0],incy,&z[0],incz);
}
    
        /// \brief Zero vector
template<class T>  void Zero(int n, Array<OneD,T> &x, const int incx)
{
    Zero(n,&x[0],incx);
}
        
        /// \brief Negate x = -x
template<class T>  void Neg( int n, Array<OneD,T> &x, const int incx)
{
    Neg(n,&x[0],incx);
}
    
template<class T> void Vlog(int n, const Array<OneD,const T> &x, const int incx, Array<OneD,T> &y, const int incy)
{
    Vlog(n, &x[0], incx, &y[0], incy);
}


template<class T> void Vexp(int n, const Array<OneD,const T> &x, const int incx, Array<OneD,T> &y, const int incy)
{
    Vexp(n, &x[0], incx, &y[0], incy);
}

template<class T> void Vpow(int n, const Array<OneD,const T> &x, const int incx, const T f, Array<OneD,T> &y, const int incy)
{
    Vpow(n, &x[0], incx, f, &y[0], incy);
}

        /// \brief sqrt y = sqrt(x)
template<class T> void Vsqrt(int n, const Array<OneD,const T> &x, const int incx, Array<OneD,T> &y, const int incy)
{
    Vsqrt(n,&x[0],incx,&y[0],incy);
}
    
        /// \brief vabs: y = |x|
template<class T> void Vabs(int n, const Array<OneD,const T> &x, const int incx, Array<OneD,T> &y, const int incy)
{
    Vabs(n,&x[0],incx,&y[0],incy);
}
    
        /********** Triad  routines  ***********************/
        
        /// \brief  vvtvp (vector times vector plus vector): z = w*x + y
template<class T> void Vvtvp(int n, const Array<OneD, const T> &w, const int incw, const Array<OneD,const T> &x, const int incx, const Array<OneD, const T> &y, const int incy, Array<OneD,T> &z, const int incz)
{
    Vvtvp(n,&w[0],incw,&x[0],incx,&y[0],incy,&z[0],incz);
}

        /// \brief  svtvp (scalar times vector plus vector): z = alpha*x + y
template<class T> void Svtvp(int n, const T alpha, const Array<OneD,const T> &x,  const int incx, const Array<OneD, const T> &y, const int incy, Array<OneD,T> &z, const int incz)
{
    Svtvp(n,alpha,&x[0],incx,&y[0],incy,&z[0],incz);
}

        /// \brief  svtvp (scalar times vector plus vector): z = alpha*x + y
template<class T> void Svtvm(int n, const T alpha, const Array<OneD,const T> &x,  const int incx, const Array<OneD, const T> &y, const int incy, Array<OneD,T> &z, const int incz)
{
    Svtvm(n,alpha,&x[0],incx,&y[0],incy,&z[0],incz);
}

        /// \brief vvtvm (vector times vector minus vector): z = w*x - y
template<class T> void Vvtvm(int n, const Array<OneD,const T> &w, const int incw, const Array<OneD,const T> &x, const int incx, const Array<OneD,const T> &y, const int incy,  Array<OneD,T> &z, const int incz)
{
    Vvtvm(n,&w[0],incw,&x[0],incx,&y[0],incy,&z[0],incz);
}
        
        /// \brief vvtvvtp (vector times vector plus vector times vector): z = v*w + y*z
template<class T> void Vvtvvtp (
            int n,
            const Array<OneD,const T> &v, int incv,
            const Array<OneD,const T> &w, int incw,
            const Array<OneD,const T> &x, int incx,
            const Array<OneD,const T> &y, int incy,
                  Array<OneD,      T> &z, int incz)
{
    Vvtvvtp(n,&v[0],incv,&w[0],incw,&x[0],incx,&y[0],incy,&z[0],incz);
}

        /// \brief svtsvtp (scalar times vector plus scalar times vector): z = alpha*x + beta*y
template<class T> void Svtsvtp(int n, const T alpha, const Array<OneD,const T> &x, const int incx, const T beta, const Array<OneD,const T> &y, const int incy,  Array<OneD,T> &z, const int incz)
{
    Svtsvtp(n,alpha,&x[0],incx,beta,&y[0],incy,&z[0],incz);
}

template<class T> void Vcopy(int n, const Array<OneD, const T> &x, int incx, Array<OneD,T> &y, int const incy)
{
    Vcopy(n,&x[0],incx,&y[0],incy);
}
    
}
#endif //VMATHARRAY_HPP
