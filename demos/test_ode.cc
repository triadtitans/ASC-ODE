#include <nonlinfunc.h>
#include <ode.h>

using namespace ASC_ode;

class RHS : public NonlinearFunction
{
  size_t DimX() const override { return 2; }
  size_t DimF() const override { return 2; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -x(0);
  }
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -1;
  }
};


int main()
{
  double tend = 1;
  double dt = 0.01;
  Vector<> x { 1, 0 };
  auto rhs = make_shared<RHS>();
  
  SolveODE_IE(tend, dt, x, rhs,
              [](double t, VectorView<double> x) { cout << "t = " << t << ", x = " << x << endl; });
}
