#ifndef Newton_h
#define Newton_h

#include <iomanip>
#include "nonlinfunc.h"
#include "lapack_interface.h"

namespace ASC_ode
{

  void NewtonSolver (std::shared_ptr<NonlinearFunction> func, VectorView<double> x,
                     double tol = 1e-10, int maxsteps = 10,
                     std::function<void(int,double,VectorView<double>)> callback = nullptr)
  {
    Vector<> res(func->DimF());
    Matrix<> fprime(func->DimF(), func->DimX());

    //std::cout << "x = " << x << std::endl;
    for (int i = 0; i < maxsteps; i++)
      {
        
        func->Evaluate(x, res);
        std::cout << "f: " << res << std::endl << std::endl;
        //std::cout<< "eval" << std::endl;
        
        func->EvaluateDeriv(x, fprime);
        break;
        //std::cout << std::setprecision(2) << "fprime = " << fprime << std::endl;
        Matrix<double> fprime_inv = inverse(fprime);

        //LapackLU<Ordering::RowMajor> lu_fprime(fprime);
        //Matrix<double> fprime_inv(lu_fprime.Inverse());
        //fprime = lu_fprime.Inverse();
        
        x -= fprime_inv*res;
        
        double err = Norm(res);
        if (callback)
          callback(i, err, x);
        if (err < tol) return;
      }

    throw std::domain_error("Newton did not converge");
  }

}

#endif
