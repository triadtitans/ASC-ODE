#include <nonlinfunc.h>
#include <ode.h>
#include <math.h>
using namespace ASC_ode;


class MassSpring : public NonlinearFunction
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

class ElectricNetwork : public NonlinearFunction
{
  size_t DimX() const override { return 2; }
  size_t DimF() const override { return 2; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = (cos(100*M_PI*x(1))- x(0));
    f(1) = 1;
  }
  
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,0) = -1;
    df(0,1) = -(100*M_PI)*sin(100*M_PI*x(1));
  }
};

int main()
{
  double tend = 32*M_PI;
  int steps = 100;
  Vector<double> y {2};
  y(0)=1;
  y(1)=0;
  auto rhs = std::make_shared<ElectricNetwork>();
  
  SolveODE_CN(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { std::cout << t << "  " << y(0) << " " << y(1) << std::endl; });


}
