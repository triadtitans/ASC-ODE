#include <nonlinfunc.h>
#include <math.h>
#include <ode.h>

using namespace ASC_ode;

// the pendulum with a length constraint

// Lagrange = -f*y + lam*(x*x+y*y-1)
// dLagrange
class dLagrangeold : public NonlinearFunction
{
  size_t DimX() const override { return 3; }
  size_t DimF() const override { return 3; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = 2*x(0)*x(2);
    f(1) = 2*x(1)*x(2) - 1;
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

class rhsLagrange : public NonlinearFunction
{
  public:
  size_t DimX() const override { return 18; }
  size_t DimF() const override { return 1; }
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    //ACHTUNG! KEINA AHNUNG OB INDEX +-1 
    VectorView<double> bview = x.Range(3,12);
    MatrixView<double> b = AsMatrix(bview,3,3);
    auto c = TransposeMatExpr(b)* b + (-1)*IdMatExpr(3);
    Vector<double> g(6);
    g(0)=c(0, 0);
    g(1)=c(1, 0);
    g(2)=c(2, 0);
    g(3)=c(1, 1);
    g(4)=c(2, 1);
    g(5)=c(2, 2);
    f(0) = x.Range(12, 18)*g;
  }
  virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
  {
    // TODO: exact differentiation
    double eps = 1e-8;
    Vector<> xl(DimX()), xr(DimX()), fl(DimF()), fr(DimF());
    for (size_t i = 0; i < DimX(); i++)
      {
        xl = x;
        xl(i) -= eps;
        xr = x;
        xr(i) += eps;
        Evaluate (xl, fl);
        Evaluate (xr, fr);
        df.Col(i) = 1/(2*eps) * (fr-fl);
      }
  }
};

class dLagrange : public NonlinearFunction
{
  size_t DimX() const override { return 18; }
  size_t DimF() const override { return 18; }
  rhsLagrange rhs;
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    Matrix<double> fmat (1, 18);
    rhs.EvaluateDeriv(x, fmat);
    f = fmat.Row(0);
  }
  virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
  {
    // TODO: exact differentiation
    double eps = 1e-8;
    Vector<> xl(DimX()), xr(DimX()), fl(DimF()), fr(DimF());
    for (size_t i = 0; i < DimX(); i++)
      {
        xl = x;
        xl(i) -= eps;
        xr = x;
        xr(i) += eps;
        Evaluate (xl, fl);
        Evaluate (xr, fr);
        df.Col(i) = 1/(2*eps) * (fr-fl);
      }
  }
};


int main()
{
  //std::cout << "main";
  // double tend = 50*2*M_PI;
  // double steps = 10000;
  // Vector<double> x { 3 };
  // x(0)=1;
  // Vector<double> dx { 3 };
  // Vector<double> ddx { 3 };
  // auto rhs = std::make_shared<dLagrange>();
  // auto mass = std::make_shared<Projector>(3, 0, 2);
  
  // SolveODE_Alpha (tend, steps, 0.8, x, dx, ddx, rhs, mass, 
  //                  // [](double t, VectorView<double> x) { cout << "t = " << t << ", x = " << x(0) << " " << x(1) << " " << x(2) << endl; }
  //                  [](double t, VectorView<double> x) { std::cout << t << " " << x(0) << " " << x(1) << " " << x(2) << std::endl; }                   
  //                  );

  double tend = 50*2*M_PI;
  double steps = 10000;
  Vector<double> x { 6 };
  x(0)=1;
  x(1)=0;
  x(2)=1;
  x(3)=-1;
  Vector<double> dx { 6 };
  Vector<double> ddx { 6 };
  auto rhs = std::make_shared<dLagrange>();
  auto mass = std::make_shared<Projector>(6, 0, 4);
  std::cout << "list:={";
  SolveODE_Alpha (tend, steps, 0.8, x, dx, ddx, rhs, mass, 
                   // [](double t, VectorView<double> x) { cout << "t = " << t << ", x = " << x(0) << " " << x(1) << " " << x(2) << endl; }
                   [](double t, VectorView<double> x) { std::cout<<std::fixed <<"{"<< t << " ,{" << x(0) << " ," << x(1) << "},{ " << x(2) << " ," << x(3) << "}}," << std::endl; }                   
                   );
  std::cout << "}";
}
