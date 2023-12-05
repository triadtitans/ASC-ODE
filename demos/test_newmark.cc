#include <nonlinfunc.h>
#include <math.h>
#include <ode.h>

using namespace ASC_ode;

class RHS : public NonlinearFunction
{
  size_t DimX() const override { return 1; }
  size_t DimF() const override { return 1; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = -x(0);
  }
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    //df(0) = -1;
  }
};


int main()
{
  double tend = 2*M_PI;
  int steps = 100;
  Vector<> x { 2 };
  x(0)=1;
  x(1)=0;
  Vector<> dx { 2 };
  dx(0)=0;
  dx(1)=0;
  auto rhs = std::make_shared<RHS>();
  auto mass = std::make_shared<IdentityFunction>(1);
  SolveODE_Newmark(tend, steps, x, dx, rhs, mass,
                   [](double t, VectorView<double> x) { std::cout << "t = " << t << ", x = " << x(0) << std::endl; }
                   );
}
