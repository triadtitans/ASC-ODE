#include <nonlinfunc.h>
#include <math.h>
#include <ode.h>

using namespace ASC_ode;

// the pendulum with a length constraint

// Lagrange = -f*y + lam*(x*x+y*y-1)
// dLagrange
class dLagrange : public NonlinearFunction
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


class dLagrangeDoublePendulum : public NonlinearFunction
{
  size_t DimX() const override { return 6; }
  size_t DimF() const override { return 6; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    //x = (x_1_x,x_1_y,x_2_x,x_2_y,l_1,l_2)
    //       0     1    2     3     4   5

    f(0) = 2*x(4)*x(0) + 2*x(5)*(x(0)-x(2));
    f(1) = -1+2*x(1)*x(4)+2*x(5)*(x(1)-x(3));
    f(2) = -2*x(5)*(x(0)-x(2));
    f(3) = -1-2*x(5)*(x(1)-x(3));
    f(4) = x(0)*x(0)+x(1)*x(1)-1;
    f(5) = (x(0)-x(2))*(x(0)-x(2))+(x(1)-x(3))*(x(1)-x(3))-1;
    
  }
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df(0,0) = 2*x(4) + 2*x(5);
    df(0,1) = 0;
    df(0,2) = -2*x(5);
    df(0,3) = 0;
    df(0,4) = 2*x(0);
    df(0,5) = 2*(x(0)-x(2));

    df(1,0) = 0;
    df(1,1) = 2*x(4)+2*x(5);
    df(1,2) = 0;
    df(1,3) = 2*x(5);
    df(1,4) = 2*x(1);
    df(1,5) = 2*(x(1)-x(3));

    df(2,0) = -2*x(5);
    df(2,1) = 0;
    df(2,2) = 2*x(5);
    df(2,3) = 0;
    df(2,4) = 0;
    df(2,5) = -2*(x(0)-x(2));

    df(3,0) = 0;
    df(3,1) = -2*x(5);
    df(3,2) = 0;
    df(3,3) = 2*x(5);
    df(3,4) = 0;
    df(3,5) = -2*(x(1)-x(3));

    df(4,0) = 2*x(0);
    df(4,1) = 2*x(1);
    df(4,2) = 0;
    df(4,3) = 0;
    df(4,4) = 0;
    df(4,5) = 0;

    df(5,0) = 2*(x(0)-x(2));
    df(5,1) = 2*(x(1)-x(3));
    df(5,2) = -2*(x(0)-x(2));
    df(5,3) = -2*(x(1)-x(3));
    df(5,4) = 0;
    df(5,5) =0;
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
  Vector<double> x (6);
  x(0)=1;
  x(1)=0;
  x(2)=1;
  x(3)=-1;
  Vector<double> dx (6);
  Vector<double> ddx (6);
  auto rhs = std::make_shared<dLagrangeDoublePendulum>();
  auto mass = std::make_shared<Projector>(6, 0, 4);
  std::cout << "list:={";
  SolveODE_Alpha (tend, steps, 0.8, x, dx, ddx, rhs, mass, 
                   // [](double t, VectorView<double> x) { cout << "t = " << t << ", x = " << x(0) << " " << x(1) << " " << x(2) << endl; }
                   [](double t, VectorView<double> x) { std::cout<<std::fixed <<"{"<< t << " ,{" << x(0) << " ," << x(1) << "},{ " << x(2) << " ," << x(3) << "}}," << std::endl; }                   
                   );
  std::cout << "}";
}
