#include <nonlinfunc.h>
#include <ode.h>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
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
  double RC=0.01;
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = (cos(100*M_PI*x(1))- x(0))*1/(RC);
    f(1) = 1;
  }
  
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,0) = -1*1/(RC);
    df(0,1) = -(100*M_PI)*sin(100*M_PI*x(1))*1/(RC);
  }
};
std::string formatVec(std::vector<double> v) { 
  std::string outstring = "";
  for( size_t i = 0; i < v.size(); ++i ) {
    outstring += std::to_string(v[i]);
    outstring += ",";
  }
  outstring.pop_back();
  return outstring;
}

int main()
{
  double tend = 0.1*M_PI;
  int steps = 100;
  Vector<double> y {2};
  y(0)=0;
  y(1)=0;
  auto rhs = std::make_shared<ElectricNetwork>();

  std::vector<double> time;
  std::vector<double> f1;
  
  SolveODE_IE(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { std::cout << t << "  " << y(0) << " " << y(1) << std::endl; });

  std::ofstream file;
  file.open("./data.txt", std::ios::trunc);

  file << "time := {" << formatVec(time) << "}" << std::endl;
  file << "function1 := {" << formatVec(f1) << "}" << std::endl;
  file << "points = Transpose[{time, function1}]\n" << "ListPlot[points, Joined -> True, PlotMarkers -> Automatic, PlotStyle -> Blue]";
  file.close();
}
