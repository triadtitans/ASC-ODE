#ifndef Newton_h
#define Newton_h

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
        //std::cout << "res ="<< res << std::endl;
        func->Evaluate(x, res);
        std::cout << "res ="<< res << std::endl;
        //std::cout << "|res| = " << Norm(res) << std::endl;
        func->EvaluateDeriv(x, fprime);
        //std::cout << "fprime = " << fprime << std::endl;
        // fprime = inverse(fprime);
        LapackLU<Ordering::RowMajor> lu_fprime(fprime);
        Matrix<double> fprime_inv(lu_fprime.Inverse());
        // fprime = lu_fprime.Inverse();
        //std::cout << "inv fprime = " << fprime << std::endl;
        x -= fprime_inv*res;
        //std::cout << "new x = " << x << std::endl;

        double err = Norm(res);
        if (callback)
          callback(i, err, x);
        if (err < tol) return;
      }

    throw std::domain_error("Newton did not converge");
  }

}

#endif
