#ifndef ODE_h
#define ODE_h

#include <functional>
#include <exception>

#include "Newton.h"


namespace ASC_ode
{

  class BlockMatVec: public NonlinearFunction {
    size_t DimX() const override { return w_ * n_; }
    size_t DimF() const override { return h_ * n_; }
    double RC = 0.01;
    int w_;
    int h_;
    int n_;
    Matrix<double> m1_;

    BlockMatVec(int h, int w, int n, MatrixView<double> m1): h_(h), w_(w), n_(n), m1_(m1) {};

    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      // A has dimension h_ x w_
      for(size_t j = 0; j < h_; j++) {
        f.Range(j * n_, (j + 1) * n_) = 0;
        for(size_t l = 0; l < w_; l++) {
            f.Range(j * n_, (j + 1) * n_) += m1_(j,l)*x.Slice(l * n_, (l + 1) * n_);
        };
      };
    }
    
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    { 
      for(size_t j = 0; j < h_; j++) {
        for(size_t l = 0; l < w_; l++) {
          df.Rows(j * h_, n_).Cols(l * w_, n_) = m1_(j,l);
        }
      }
    }
  };
  
  // implicit Euler method for dy/dt = rhs(y)
  void SolveODE_IE(double tend, int steps,
                   VectorView<double> y, std::shared_ptr<NonlinearFunction> rhs,
                   std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    auto yold = std::make_shared<ConstantFunction>(y);
    auto ynew = std::make_shared<IdentityFunction>(y.Size());
    auto equ = ynew-yold - dt * rhs;

    double t = 0;
    for (int i = 0; i < steps; i++)
      {
        NewtonSolver (equ, y);
        yold->Set(y);
        t += dt;
        if (callback) callback(t, y);
      }
  }

  void SolveODE_EE(double tend, int steps,
                   VectorView<double> y, std::shared_ptr<NonlinearFunction> rhs,
                   std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    Vector<double> dy(rhs->DimF());
    double t = 0;
    for (int i = 0; i < steps; i++)
      {
        rhs->Evaluate(y,dy);
        y += dt*dy;
        t += dt;
        if (callback) callback(t, y);
      }
  }

 void SolveODE_CN(double tend, int steps,
                   VectorView<double> y, std::shared_ptr<NonlinearFunction> rhs,
                   std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    auto yold = std::make_shared<ConstantFunction>(y);
    Vector<double> fold(rhs->DimF());
    rhs->Evaluate(y,fold);
    auto rhsold = std::make_shared<ConstantFunction>(fold);
    auto ynew = std::make_shared<IdentityFunction>(y.Size());
    auto equ = ynew-yold - dt * 1/2 * (rhsold + rhs);

    double t = 0;
    for (int i = 0; i < steps; i++)
      {
        NewtonSolver (equ, y);
        yold->Set(y);
        rhs->Evaluate(y,fold);
        rhsold->Set(fold);
        
        t += dt;
        if (callback) callback(t, y);
      }
  }
  

  
  
  
  // Newmark and generalized alpha:
  // https://miaodi.github.io/finite%20element%20method/newmark-generalized/
  
  // Newmark method for  mass*d^2x/dt^2 = rhs
  void SolveODE_Newmark(double tend, int steps,
                        VectorView<double> x, VectorView<double> dx,
                        std::shared_ptr<NonlinearFunction> rhs,   
                        std::shared_ptr<NonlinearFunction> mass,  
                        std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    double gamma = 0.5;
    double beta = 0.25;

    Vector<> a(x.Size());
    Vector<> v(x.Size());

    auto xold = std::make_shared<ConstantFunction>(x);
    auto vold = std::make_shared<ConstantFunction>(dx);
    auto aold = std::make_shared<ConstantFunction>(x);
    rhs->Evaluate (xold->Get(), aold->Get());
    
    auto anew = std::make_shared<IdentityFunction>(a.Size());
    auto vnew = vold + dt*((1-gamma)*aold+gamma*anew);
    auto xnew = xold + dt*vold + dt*dt/2 * ((1-2*beta)*aold+2*beta*anew);    

    auto equ = Compose(mass, anew) - Compose(rhs, xnew);

    double t = 0;
    for (int i = 0; i < steps; i++)            
      {
        NewtonSolver (equ, a);
        xnew -> Evaluate (a, x);
        vnew -> Evaluate (a, v);

        xold->Set(x);
        vold->Set(v);
        aold->Set(a);
        t += dt;
        if (callback) callback(t, x);
      }
    dx = v;
  }




  // Generalized alpha method for M d^2x/dt^2 = rhs
  void SolveODE_Alpha (double tend, int steps, double rhoinf,
                       VectorView<double> x, VectorView<double> dx, VectorView<double> ddx,
                       std::shared_ptr<NonlinearFunction> rhs,   
                       std::shared_ptr<NonlinearFunction> mass,  
                       std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double dt = tend/steps;
    double alpham = (2*rhoinf-1)/(rhoinf+1);
    double alphaf = rhoinf/(rhoinf+1);
    double gamma = 0.5-alpham+alphaf;
    double beta = 0.25 * (1-alpham+alphaf)*(1-alpham+alphaf);

    Vector<> a(x.Size());
    Vector<> v(x.Size());

    auto xold = std::make_shared<ConstantFunction>(x);
    auto vold = std::make_shared<ConstantFunction>(dx);
    auto aold = std::make_shared<ConstantFunction>(ddx);
    // rhs->Evaluate (xold->Get(), aold->Get()); // solve with M ???
    
    auto anew = std::make_shared<IdentityFunction>(a.Size());
    auto vnew = vold + dt*((1-gamma)*aold+gamma*anew);
    auto xnew = xold + dt*vold + dt*dt/2 * ((1-2*beta)*aold+2*beta*anew);    

    // auto equ = Compose(mass, (1-alpham)*anew+alpham*aold) - Compose(rhs, (1-alphaf)*xnew+alphaf*xold);
    auto equ = Compose(mass, (1-alpham)*anew+alpham*aold) - (1-alphaf)*Compose(rhs,xnew) - alphaf*Compose(rhs, xold);

    double t = 0;
    a = ddx;

    for (int i = 0; i < steps; i++)
      {
        NewtonSolver (equ, a);
        xnew -> Evaluate (a, x);
        vnew -> Evaluate (a, v);

        xold->Set(x);
        vold->Set(v);
        aold->Set(a);
        t += dt;
        if (callback) callback(t, x);
      }
    dx = v;
    ddx = a;
  }

  

}


#endif
