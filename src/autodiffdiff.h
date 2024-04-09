#ifndef FILE_AUTODIFFDIFF
#define FILE_AUTODIFFDIFF
#include <iostream>
#include <math.h>

using namespace std;

namespace ASC_ode
{
    template <size_t D, typename SCAL = double>
    class AutoDiffDiff
    {
        SCAL val;
        SCAL dval[D?D:1] ={0};
        SCAL ddval[D?D*D:1]={0};
    public:

        typedef AutoDiffDiff<D, SCAL> TELEM;


        /// elements are undefined
        AutoDiffDiff  () throw() { ; }

        /// copy constructor
        AutoDiffDiff  (const AutoDiffDiff & ad2) throw()
        {
            val = ad2.val;
            for (size_t i = 0; i < D; i++)
            dval[i] = ad2.dval[i];
            for (size_t i = 0; i < D*D; i++)
            ddval[i] = ad2.ddval[i];
        }

        /// initial object with constant value
        AutoDiffDiff  (SCAL aval) throw()
        {
            val = aval;
            for (size_t i = 0; i < D; i++)
            dval[i] = 0;
            for (size_t i = 0; i < D*D; i++)
            ddval[i] = 0;
        }

        /// initial object with value and derivative
        /*
        AutoDiffDiff  (const AutoDiffVec<D, SCAL> & ad2) throw()
        {
            val = ad2.Val();
            for (size_t i = 0; i < D; i++)
            dval[i] = ad2.DValue(i);
            for (size_t i = 0; i < D*D; i++)
            ddval[i] = 0;
        }
        */

        /// init object with (val, e_diffindex)
        AutoDiffDiff  (SCAL aval, size_t diffindex)  throw()
        {
            val = aval;
            for (size_t i = 0; i < D; i++)
            dval[i] = 0;
            for (size_t i = 0; i < D*D; i++)
            ddval[i] = 0;
            dval[diffindex] = 1;
        }

        AutoDiffDiff (SCAL aval, const SCAL * grad)
        {
            val = aval;
            LoadGradient (grad);
            for (size_t i = 0; i < D*D; i++)
            ddval[i] = 0;
        }

        AutoDiffDiff (SCAL aval, const SCAL * grad, const SCAL * hesse)
        {
            val = aval;
            LoadGradient (grad);
            LoadHessian (hesse);
        }

        /// assign constant value
        AutoDiffDiff & operator= (SCAL aval) throw()
        {
            val = aval;
            for (size_t i = 0; i < D; i++)
                dval[i] = 0;
            for (size_t i = 0; i < D*D; i++)
                ddval[i] = 0;
            return *this;
        }

        void StoreGradient (SCAL * p) const 
        {
            for (size_t i = 0; i < D; i++)
                p[i] = dval[i];
        }

        void LoadGradient (const SCAL * p) 
        {
            for (size_t i = 0; i < D; i++)
                dval[i] = p[i];
        }

        void StoreHessian (SCAL * p) const 
        {
            for (size_t i = 0; i < D*D; i++)
                p[i] = ddval[i];
        }

        void LoadHessian (const SCAL * p) 
        {
            for (size_t i = 0; i < D*D; i++)
                ddval[i] = p[i];
        }

        /// returns value
        SCAL Value() const throw() { return val; }

        /// returns partial derivative
        SCAL DValue (size_t i) const throw() { return dval[i]; }

        /*
        AutoDiffVec<D,SCAL> DValueAD (size_t i) const
        {
            AutoDiffVec<D,SCAL> r(dval[i]);
            for (size_t j = 0; j < D; j++)
            r.DVal(j) = ddval[i*D+j];
            return r;
        }
        */
        
        /// returns partial derivative
        SCAL DDValue (size_t i) const throw() { return ddval[i]; }

        /// returns partial derivative
        SCAL DDValue (size_t i, size_t j) const throw() { return ddval[i*D+j]; }

        /// access value
        SCAL & Value() throw() { return val; }

        /// accesses partial derivative 
        SCAL & DValue (size_t i) throw() { return dval[i]; }

        /// accesses partial derivative 
        SCAL & DDValue (size_t i) throw() { return ddval[i]; }

        /// accesses partial derivative 
        SCAL & DDValue (size_t i, size_t j) throw() { return ddval[i*D+j]; }

        /*
        explicit operator AutoDiffVec<D,SCAL> () const
        { return AutoDiffVec<D,SCAL> (val, &dval[0]); }
        */

        /// add autodiffdiff object
        AutoDiffDiff<D, SCAL> & operator+= (const AutoDiffDiff<D, SCAL> & y) throw()
        {
            val += y.val;
            for (size_t i = 0; i < D; i++)
                dval[i] += y.dval[i];
            for (size_t i = 0; i < D*D; i++)
                ddval[i] += y.ddval[i];
            return *this;
        }

        /// subtract autodiffdiff object
        AutoDiffDiff<D, SCAL> & operator-= (const AutoDiffDiff<D, SCAL> & y) throw()
        {
            val -= y.val;
            for (size_t i = 0; i < D; i++)
                dval[i] -= y.dval[i];
            for (size_t i = 0; i < D*D; i++)
                ddval[i] -= y.ddval[i];
            return *this;
        }

        /// multiply with autodiffdiff object
        AutoDiffDiff<D, SCAL> & operator*= (const AutoDiffDiff<D, SCAL> & y) throw()
        {
            for (size_t i = 0; i < D*D; i++) {
                ddval[i] = val * y.ddval[i] + y.val * ddval[i];
            }

            for (size_t i = 0; i < D; i++){
                for (size_t j = 0; j < D; j++){
                ddval[i*D+j] += dval[i] * y.dval[j] + dval[j] * y.dval[i];
                }
            }

            for (size_t i = 0; i < D; i++)
            {
                dval[i] *= y.val;
                dval[i] += val * y.dval[i];
            }
            val *= y.val;
            return *this;
        }

        /// multiply with scalar
        AutoDiffDiff<D, SCAL> & operator*= (const SCAL & y) throw()
        {
            for ( size_t i = 0; i < D*D; i++ ) {
                ddval[i] *= y;
            }
            for (size_t i = 0; i < D; i++) {
                dval[i] *= y;
            }
            val *= y;
            return *this;
        }

        /// divide by scalar
        AutoDiffDiff<D, SCAL> & operator/= (const SCAL & y) throw()
        {
            SCAL iy = 1.0 / y;
            for ( size_t i = 0; i < D*D; i++ ) {
                ddval[i] *= iy;
            }
            for (size_t i = 0; i < D; i++) {
                dval[i] *= iy;
            }
            val *= iy;
            return *this;
        }
    };

    template<size_t D, typename SCAL, typename SCAL2>
    inline AutoDiffDiff<D, SCAL> operator+ (SCAL2 x, const AutoDiffDiff<D, SCAL> & y) throw()
    {
        AutoDiffDiff<D, SCAL> res;
        res.Value() = x+y.Value();
        for (size_t i = 0; i < D; i++) {
            res.DValue(i) = y.DValue(i);
        }
        for (size_t i = 0; i < D*D; i++) {
            res.DDValue(i) = y.DDValue(i);
        }
        return res;
    }

    template<size_t D, typename SCAL>
    inline AutoDiffDiff<D, SCAL> operator+ (const AutoDiffDiff<D, SCAL> & x, const AutoDiffDiff<D, SCAL> & y) throw()
    {
        AutoDiffDiff<D, SCAL> res;
        res.Value() = x.Value()+y.Value();
        for (size_t i = 0; i < D; i++) {
            res.DValue(i) = x.DValue(i) + y.DValue(i);
        }
        for (size_t i = 0; i < D*D; i++) {
            res.DDValue(i) = x.DDValue(i) + y.DDValue(i);
        }
        return res;
    }

    ///
    template<size_t D, typename SCAL, typename SCAL2>
    inline AutoDiffDiff<D, SCAL> operator+ (const AutoDiffDiff<D, SCAL> & y, SCAL2 x) throw()
    {
        AutoDiffDiff<D, SCAL> res;
        res.Value() = x+y.Value();
        for (size_t i = 0; i < D; i++) {
            res.DValue(i) = y.DValue(i);
        }
        for (size_t i = 0; i < D*D; i++) {
            res.DDValue(i) = y.DDValue(i);
        }
        return res;
    }

    template<size_t D, typename SCAL>
    inline AutoDiffDiff<D, SCAL> operator- (const AutoDiffDiff<D, SCAL> & x, const AutoDiffDiff<D, SCAL> & y) throw()
    {
        AutoDiffDiff<D, SCAL> res;
        res.Value() = x.Value()-y.Value();
        for (size_t i = 0; i < D; i++) {
            res.DValue(i) = x.DValue(i) - y.DValue(i);
        }
        for (size_t i = 0; i < D*D; i++) {
            res.DDValue(i) = x.DDValue(i) - y.DDValue(i);
        }
        return res;
    }


    ///
    template<size_t D, typename SCAL>
    inline AutoDiffDiff<D, SCAL> operator- (const AutoDiffDiff<D, SCAL> & x) throw()
    {
        AutoDiffDiff<D, SCAL> res;
        res.Value() = -x.Value();
        for (size_t i = 0; i < D; i++) {
            res.DValue(i) = -x.DValue(i);
        }
        for (size_t i = 0; i < D*D; i++) {
            res.DDValue(i) = -x.DDValue(i);
        }
        return res;
    }

    ///
    template<size_t D, typename SCAL, typename SCAL2>
    inline AutoDiffDiff<D, SCAL> operator- (const AutoDiffDiff<D, SCAL> & x, SCAL2 y) throw()
    {
        AutoDiffDiff<D, SCAL> res;
        res.Value() = x.Value()-y;
        for (size_t i = 0; i < D; i++) {
            res.DValue(i) = x.DValue(i);
        }
        for (size_t i = 0; i < D*D; i++) {
            res.DDValue(i) = x.DDValue(i);
        }
        return res;
    }

    ///
    template<size_t D, typename SCAL, typename SCAL2>
    inline AutoDiffDiff<D, SCAL> operator- (SCAL2 x, const AutoDiffDiff<D, SCAL> & y) throw()
    {
        AutoDiffDiff<D, SCAL> res;
        res.Value() = x-y.Value();
        for (size_t i = 0; i < D; i++) {
            res.DValue(i) = -y.DValue(i);
        }
        for (size_t i = 0; i < D*D; i++) {
            res.DDValue(i) = -y.DDValue(i);
        }
        return res;
    }


    ///
    template<size_t D, typename SCAL, typename SCAL2>
    inline AutoDiffDiff<D, SCAL> operator* (SCAL2 x, const AutoDiffDiff<D, SCAL> & y) throw()
    {
        AutoDiffDiff<D, SCAL> res;
        res.Value() = x*y.Value();
        for (size_t i = 0; i < D; i++) {
            res.DValue(i) = x*y.DValue(i);
        }
        for (size_t i = 0; i < D*D; i++) {
            res.DDValue(i) = x*y.DDValue(i);
        }
        return res;
    }

    ///
    template<size_t D, typename SCAL, typename SCAL2>
    inline AutoDiffDiff<D, SCAL> operator* (const AutoDiffDiff<D, SCAL> & y, SCAL2 x) throw()
    {
        AutoDiffDiff<D, SCAL> res;
        res.Value() = x*y.Value();
        for (size_t i = 0; i < D; i++) {
            res.DValue(i) = x*y.DValue(i);
        }
        for (size_t i = 0; i < D*D; i++) {
            res.DDValue(i) = x*y.DDValue(i);
        }
        return res;
    }

    ///
    template<size_t D, typename SCAL>
    inline AutoDiffDiff<D, SCAL> operator* (const AutoDiffDiff<D, SCAL> & x, const AutoDiffDiff<D, SCAL> & y) throw()
    {
        AutoDiffDiff<D, SCAL> res;
        SCAL hx = x.Value();
        SCAL hy = y.Value();

        res.Value() = hx*hy;
        for (size_t i = 0; i < D; i++) {
            res.DValue(i) = hx*y.DValue(i) + hy*x.DValue(i);
        }

        for (size_t i = 0; i < D; i++) {
            for (size_t j = 0; j < D; j++) {
            res.DDValue(i,j) = hx * y.DDValue(i,j) + hy * x.DDValue(i,j)
            + x.DValue(i) * y.DValue(j) + x.DValue(j) * y.DValue(i);
            }
        }

        return res;
    }



    template<size_t D, typename SCAL>
    inline AutoDiffDiff<D, SCAL> Inv (const AutoDiffDiff<D, SCAL> & x)
    {
        AutoDiffDiff<D, SCAL> res(1.0 / x.Value());
        for (size_t i = 0; i < D; i++) {
            res.DValue(i) = -x.DValue(i) / (x.Value() * x.Value());
        }

        SCAL fac1 = 2/(x.Value()*x.Value()*x.Value());
        SCAL fac2 = 1/sqr(x.Value());
        for (size_t i = 0; i < D; i++) {
            for (size_t j = 0; j < D; j++) {
            res.DDValue(i,j) = fac1*x.DValue(i)*x.DValue(j) - fac2*x.DDValue(i,j);
            }
        }
        return res;
    }


    template<size_t D, typename SCAL>
    inline AutoDiffDiff<D, SCAL> operator/ (const AutoDiffDiff<D, SCAL> & x, const AutoDiffDiff<D, SCAL> & y)
    {
        return x * Inv (y);
    }

    template<size_t D, typename SCAL, typename SCAL2>
    inline AutoDiffDiff<D, SCAL> operator/ (const AutoDiffDiff<D, SCAL> & x, SCAL2 y)
    {
        return (1/y) * x;
    }

    template<size_t D, typename SCAL, typename SCAL2>
    inline AutoDiffDiff<D, SCAL> operator/ (SCAL2 x, const AutoDiffDiff<D, SCAL> & y)
    {
        return x * Inv(y);
    }


    template<size_t D, typename SCAL>
    inline AutoDiffDiff<D, SCAL> sqrt (const AutoDiffDiff<D, SCAL> & x)
    {
        AutoDiffDiff<D, SCAL> res;
        res.Value() = std::sqrt(x.Value());
        for (size_t j = 0; j < D; j++) {
            if (x.DValue(j) == 0) {
                res.DValue(j) = SCAL{0.0};
            } else {
                res.DValue(j) = 0.5 / res.Value() * x.DValue(j);
            }
        }

        
        for (size_t i = 0; i < D; i++) {
            for (size_t j = 0; j < D; j++) {
                if((x.DDValue(i,j)+x.DValue(i)*x.DValue(j)) == 0) {
                    res.DDValue(i, j) = SCAL{0.};
                } else {
                    0.5/res.Value() * x.DDValue(i,j) - 0.25 / (x.Value()*res.Value()) * x.DValue(i) * x.DValue(j);
                }
            }
        }

        return res;
    }

    // df(u)/dx  = exp(x) * du/dx
    // d^2 f(u) / dx^2 = exp(x) * (du/dx)^2 + exp(x) * d^2u /dx^2
    template <size_t D, typename SCAL>
    inline AutoDiffDiff<D, SCAL> exp (AutoDiffDiff<D, SCAL> x)
    {
        AutoDiffDiff<D, SCAL> res;
        res.Value() = std::exp(x.Value());
        for (size_t k = 0; k < D; k++) {
            res.DValue(k) = x.DValue(k) * res.Value();
        }
        for (size_t k = 0; k < D; k++) {
            for (size_t l = 0; l < D; l++) {
            res.DDValue(k,l) = (x.DValue(k) * x.DValue(l)+x.DDValue(k,l)) * res.Value();
            }
        }
        return res;
    }

    template <size_t D, typename SCAL>
    inline AutoDiffDiff<D, SCAL> log (AutoDiffDiff<D, SCAL> x)
    {
        AutoDiffDiff<D, SCAL> res;
        res.Value() = std::log(x.Value());
        SCAL xinv = 1.0/x.Value();
        for (size_t k = 0; k < D; k++) {
            res.DValue(k) = x.DValue(k) * xinv;
        }
        for (size_t k = 0; k < D; k++) {
            for (size_t l = 0; l < D; l++) {
            res.DDValue(k,l) = -xinv*xinv*x.DValue(k) * x.DValue(l) + xinv * x.DDValue(k,l);
            }
        }
        return res;
    }

    using std::pow;
    template <size_t D, typename SCAL>
    inline AutoDiffDiff<D,SCAL> pow (AutoDiffDiff<D,SCAL> x, AutoDiffDiff<D,SCAL> y )
    {
        return exp(log(x)*y);
    }

    using std::pow;
    template <size_t D, typename SCAL>
    inline AutoDiffDiff<D, SCAL> pow (AutoDiffDiff<D, SCAL> x, SCAL y)
    {
        AutoDiffDiff<D, SCAL> y_diff(y);
        
        return pow(x, y_diff);
    }

    template <size_t D, typename SCAL>
    inline AutoDiffDiff<D, SCAL> sin (AutoDiffDiff<D, SCAL> x)
    {
        AutoDiffDiff<D, SCAL> res;
        SCAL s = sin(x.Value());
        SCAL c = cos(x.Value());
        
        res.Value() = s;
        for (size_t k = 0; k < D; k++) {
            res.DValue(k) = x.DValue(k) * c;
        }
        for (size_t k = 0; k < D; k++) {
            for (size_t l = 0; l < D; l++) {
            res.DDValue(k,l) = -s * x.DValue(k) * x.DValue(l) + c * x.DDValue(k,l);
            }
        }
        return res;
    }


    template <size_t D, typename SCAL>
    inline AutoDiffDiff<D, SCAL> cos (AutoDiffDiff<D, SCAL> x)
    {
        AutoDiffDiff<D, SCAL> res;
        SCAL s = sin(x.Value());
        SCAL c = cos(x.Value());
        
        res.Value() = c;
        for (size_t k = 0; k < D; k++) {
            res.DValue(k) = -s * x.DValue(k);
        }
        for (size_t k = 0; k < D; k++) {
            for (size_t l = 0; l < D; l++) {
            res.DDValue(k,l) = -c * x.DValue(k) * x.DValue(l) - s * x.DDValue(k,l);
            }
        }
        return res;
    }

    template <size_t D, typename SCAL>
    inline AutoDiffDiff<D, SCAL> tan (AutoDiffDiff<D, SCAL> x)
    { return sin(x) / cos(x); }
}

#endif