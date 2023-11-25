#include <nonlinfunc.h>
#include <ode.h>

using namespace ASC_ode;


// Lagrange = -f*y + lam*(x*x+y*y-1)
// dLagrange
class dLagrangef : public NonlinearFunction
{
  size_t DimX() const override { return 3; }
  size_t DimF() const override { return 3; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = 0;
    f(1) = -1;
    f(2) = 0;
  }
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
  }
};

class dLagrangeg : public NonlinearFunction
{
  size_t DimX() const override { return 3; }
  size_t DimF() const override { return 3; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = 2*x(0)*x(2);
    f(1) = 2*x(1)*x(2);
    f(2) = x(0)*x(0)+x(1)*x(1)-1;
  }
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df(0,0) = 2*x(2);
    df(0,1) = 0;
    df(0,2) = 2*x(0);

    df(1,0) = 0;
    df(1,1) = 2*x(2);
    df(1,2) = 2*x(1);

    df(2,0) = 2*x(0);
    df(2,1) = 2*x(1);
    df(2,2) = 0;
  }
};




int main()
{
  double tend = 2*M_PI;
  double dt = tend/1000;
  Vector<> x { 1, 0, 0 };
  
  SolveODE_Shake(tend, dt, x, 2, make_shared<dLagrangef>(), make_shared<dLagrangeg>(), 
                 [](double t, VectorView<double> x) { cout << "t = " << t
                                                           << ", x = " << x(0) << " " << x(1) << " " << x(2)
                                                           << ", |x,y| = " << x(0)*x(0)+x(1)*x(1) 
                                                           << endl; });
}
