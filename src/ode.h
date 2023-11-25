#ifndef ODE_h
#define ODE_h

#include <functional>
#include <exception>
#include <calcinverse.hpp>


namespace ASC_ode
{

  void NewtonSolver (shared_ptr<NonlinearFunction> func, VectorView<double> x,
                     double tol = 1e-10, int maxsteps = 10,
                     std::function<void(int,double,VectorView<double>)> callback = nullptr)
  {
    Vector<> res(func->DimF());
    Matrix<> fprime(func->DimF(), func->DimX());

    for (int i = 0; i < maxsteps; i++)
      {
        func->Evaluate(x, res);
        // cout << "|res| = " << L2Norm(res) << endl;
        func->EvaluateDeriv(x, fprime);
        CalcInverse(fprime);
        x -= fprime*res;

        double err= L2Norm(res);
        if (callback)
          callback(i, err, x);
        if (err < tol) return;
      }

    throw std::domain_error("Newton did not converge");
  }


  // implicit Euler method for dx/dt = rhs(x)
  void SolveODE_IE(double tend, double dt,
                   VectorView<double> x, shared_ptr<NonlinearFunction> rhs,
                   std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    auto xold = make_shared<ConstantFunction>(x);
    auto xnew = make_shared<IdenticFunction>(x.Size());
    auto equ = xnew-xold - dt * rhs;

    double t = 0;
    while (t < tend)
      {
        NewtonSolver (equ, x);
        xold->Set(x);
        if (callback) callback(t, x);
        t += dt;
      }
  }

  // Vertlet method for d^2x/dt^2 = rhs
  void SolveODE_Verlet(double tend, double dt,
                       VectorView<double> x, VectorView<double> dx,
                       shared_ptr<NonlinearFunction> rhs,   // x->f(x)
                       std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    auto xold = make_shared<ConstantFunction>(x);
    auto xoldold = make_shared<ConstantFunction>(Vector<double>(x-dt*dx));    
    auto xnew = make_shared<IdenticFunction>(x.Size());
    auto rhsold = make_shared<ConstantFunction>(x);        
    auto equ = xnew-2*xold+xoldold - dt*dt * rhsold;
    
    double t = 0;
    while (t < tend)
      {
        rhs->Evaluate(xold->Get(), rhsold->Get());

        NewtonSolver (equ, x);
        xoldold->Set(xold->Get());
        xold->Set(x);
        if (callback) callback(t, x);
        t += dt;
      }
    dx = 1/dt * (xold->Get()-xoldold->Get());
  }

  // Shake algorithm = Vertlet+implicit constaints
  // d^2x/dt`2 + (dg/dt)^T lam = f
  // with g(x(t))) = 0
  // have to provide function dg : (x, lambda) -> dg/dx, dg/dlam
  void SolveODE_Shake(double tend, double dt,
                      VectorView<double> xlam, size_t nx,
                      shared_ptr<NonlinearFunction> dLagrangef,
                      shared_ptr<NonlinearFunction> dLagrangeg,                      
                      std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    auto xnew = make_shared<IdenticFunction>(xlam.Size());
    auto xold = make_shared<ConstantFunction>(xlam);
    auto xoldold = make_shared<ConstantFunction>(xlam);    
    auto dLagrangef_old = make_shared<ConstantFunction>(xlam);        
    
    // auto ddx = make_shared<ProjectFunction> (xnew-2*xold+xoldold, 0, nx);
    auto ddx = Compose(make_shared<Projector>(xlam.Size(), 0, nx), xnew-2*xold+xoldold);
    auto equ = ddx - dt*dt*(dLagrangeg+dLagrangef_old);
    double t = 0;
    while (t < tend)
      {
        dLagrangef->Evaluate(xold->Get(), dLagrangef_old->Get());        
        NewtonSolver (equ, xlam);
        xoldold->Set(xold->Get());        
        xold->Set(xlam);
        if (callback) callback(t, xlam);
        t += dt;
      }
  }
  

  
  // Newmark and generalized alpha:
  // https://miaodi.github.io/finite%20element%20method/newmark-generalized/
  

  // Newmark method for d^2x/dt^2 = rhs
  void SolveODE_Newmark(double tend, double dt,
                        VectorView<double> x, VectorView<double> dx,
                        shared_ptr<NonlinearFunction> rhs,   
                        shared_ptr<NonlinearFunction> mass,  
                        std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double gamma = 0.5;
    double beta = 0.25;

    Vector<> a(x.Size());
    Vector<> v(x.Size());

    auto xold = make_shared<ConstantFunction>(x);
    auto vold = make_shared<ConstantFunction>(dx);
    auto aold = make_shared<ConstantFunction>(x);
    rhs->Evaluate (xold->Get(), aold->Get());
    
    auto anew = make_shared<IdenticFunction>(a.Size());
    auto vnew = vold + dt*((1-gamma)*aold+gamma*anew);
    auto xnew = xold + dt*vold + dt*dt/2 * ((1-2*beta)*aold+2*beta*anew);    

    auto equ = Compose(mass, anew) - Compose(rhs, xnew);

    double t = 0;
    while (t < tend)
      {
        NewtonSolver (equ, a);
        xnew -> Evaluate (a, x);
        vnew -> Evaluate (a, v);

        xold->Set(x);
        vold->Set(v);
        aold->Set(a);
        if (callback) callback(t, x);
        t += dt;
      }
    dx = v;
  }




  // Generalized alpha method for M d^2x/dt^2 = rhs
  void SolveODE_Alpha (double tend, double dt, double rhoinf,
                       VectorView<double> x, VectorView<double> dx, VectorView<double> ddx,
                       shared_ptr<NonlinearFunction> rhs,   
                       shared_ptr<NonlinearFunction> mass,  
                       std::function<void(double,VectorView<double>)> callback = nullptr)
  {
    double alpham = (2*rhoinf-1)/(rhoinf+1);
    double alphaf = rhoinf/(rhoinf+1);
    double gamma = 0.5-alpham+alphaf;
    double beta = 0.25 * (1-alpham+alphaf)*(1-alpham+alphaf);

    Vector<> a(x.Size());
    Vector<> v(x.Size());

    auto xold = make_shared<ConstantFunction>(x);
    auto vold = make_shared<ConstantFunction>(dx);
    auto aold = make_shared<ConstantFunction>(ddx);
    // rhs->Evaluate (xold->Get(), aold->Get()); // solve with M ???
    
    auto anew = make_shared<IdenticFunction>(a.Size());
    auto vnew = vold + dt*((1-gamma)*aold+gamma*anew);
    auto xnew = xold + dt*vold + dt*dt/2 * ((1-2*beta)*aold+2*beta*anew);    

    // auto equ = Compose(mass, (1-alpham)*anew+alpham*aold) - Compose(rhs, (1-alphaf)*xnew+alphaf*xold);
    auto equ = Compose(mass, (1-alpham)*anew+alpham*aold) - (1-alphaf)*Compose(rhs,xnew) - alphaf*Compose(rhs, xold);

    double t = 0;
    a = ddx;
    while (t < tend)
      {
        NewtonSolver (equ, a);
        xnew -> Evaluate (a, x);
        vnew -> Evaluate (a, v);

        xold->Set(x);
        vold->Set(v);
        aold->Set(a);
        if (callback) callback(t, x);
        t += dt;
      }
    dx = v;
    ddx = a;
  }

  

}


#endif
