#include "../src/autodiffdiff.h"
#include <iostream>
#include <math.h>
#include "../src/nonlinfunc.h"
#include "matrix.h"

constexpr size_t dim_per_body = 18;



using namespace std;

using namespace ASC_ode;

class RhsRigidBody : public NonlinearFunction
{
  public:
  //Just a reminder
  RhsRigidBody(){}; //RigidBody b
  size_t DimX() const  { return 18; }
  size_t DimF() const  { return 1; }

  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    
    Vector<AutoDiffDiff<18, double>> x_diff(this->DimX());

    Vector<AutoDiffDiff<18, double>> f_diff(this->DimF());
    
    for (size_t j = 0; j < this->DimX(); j++) 
    {
      x_diff(j).Value() = x(j);
      x_diff(j).DValue(j) = 1;
    }

    Matrix<AutoDiffDiff<dim_per_body, double>> b (3, 3);
    // extract B from Q
    b.Row(0) = x_diff.Range(1, 4);
    b.Row(1) = x_diff.Range(5, 8);
    b.Row(2) = x_diff.Range(9, 12);
    // b must be orthonormal
    Matrix<double> eye = Diagonal(3, 1);
    auto c = TransposeMatExpr(b)* b - eye;
    Vector<AutoDiffDiff<dim_per_body, double>> g(6);
    g(0)=c(0, 0);
    g(1)=c(1, 0);
    g(2)=c(2, 0);
    g(3)=c(1, 1);
    g(4)=c(2, 1);
    g(5)=c(2, 2);
    f_diff(0) = x.Range(12, 18)*g;

    for (size_t j = 0; j < this->DimX(); j++) {
      f(j) = f_diff(0).DValue(j);
    }
    }

  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    Vector<AutoDiffDiff<18, double>> x_diff(this->DimX());

    Vector<AutoDiffDiff<18, double>> f_diff(this->DimF());
    
    for (size_t j = 0; j < this->DimX(); j++) 
    {
      x_diff(j).Value() = x(j);
      x_diff(j).DValue(j) = 1;
    }

    Matrix<AutoDiffDiff<dim_per_body, double>> b (3, 3);
    // extract B from Q
    b.Row(0) = x_diff.Range(1, 4);
    b.Row(1) = x_diff.Range(5, 8);
    b.Row(2) = x_diff.Range(9, 12);
    // b must be orthonormal
    Matrix<double> eye = Diagonal(3, 1);
    auto c = Transpose(b)* b - eye;
    Vector<AutoDiffDiff<dim_per_body, double>> g(6);
    g(0)=c(0, 0);
    g(1)=c(1, 0);
    g(2)=c(2, 0);
    g(3)=c(1, 1);
    g(4)=c(2, 1);
    g(5)=c(2, 2);
    f_diff(0) = x_diff.Range(12, 18)*g;

    for (size_t j = 0; j < this->DimX(); j++) {
      for (size_t k = 0; k < this->DimX(); k++) {
          df(j, k) = f_diff(0).DDValue(j, k);
      } 
    }
  }
};

template<typename T>
T f(T x) 
{
    return sin(sqrt(exp(x)+7)/2);
}


int main() 
{
    
    RhsRigidBody body;

    Vector<double> x(18);
    x = 1;
    Vector<double> df(18);
    Vector<double> df1(18);
    Matrix<double> ddf(18, 18);
    Matrix<double> ddf1(18, 18);

    std::shared_ptr<RhsRigidBody> rhs = std::make_shared<RhsRigidBody>();
    std::shared_ptr<Derivative> dlagrange =std::make_shared<Derivative>(rhs);
    

    dlagrange->Evaluate(x, df);
    //dlagrange->Evaluate(x, df1);


    dlagrange->EvaluateDeriv(x, ddf);
    //dlagrange->EvaluateDeriv(x, ddf1);
    //body.EvaluatDeriv1(x, df1);
    //body.EvaluateDerivDeriv1(x, ddf1);

    //cout << df1 << "\n";
    cout << df;
    //cout << df1;
    cout << "hessian n" << endl;
    //cout << ddf1 << "\n";
    cout << ddf <<"\n";
    //cout << ddf1;
    
    
}


