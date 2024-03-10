#include <nonlinfunc.h>
#include <math.h>
#include <matrix.h>
#include <vector.h>
#include <ode.h>

class RhsRigidBody;

using namespace ASC_bla;
using namespace ASC_ode;

class RigidBody {
  Vector<double> q_;
  Vector<double> dq_;
  Vector<double> ddq_;
  Matrix<double> mass_;
  int dim_;
  std::shared_ptr<LinearFunction> mass_function;
public:
  template <typename T>
  RigidBody(MatrixExpr<T>& m,Vector<double> q,Vector<double> dq,Vector<double> ddq)
        : mass_(m), dim_(m.Height()), mass_function(std::make_shared<LinearFunction>(mass_)),
          q_(q),dq_(dq),ddq_(ddq){
    if(m.Width() != dim_) throw std::invalid_argument("Mass matrix must be square");
    if(q.Size() != dim_) throw std::invalid_argument("q Vector must match mass matrix");
    if(dq.Size() != dim_) throw std::invalid_argument("q Vector must match mass matrix");
    if(ddq.Size() != dim_) throw std::invalid_argument("q Vector must match mass matrix");
  }
  void simulate(double tend, double steps, std::function<void(double,VectorView<double>)> callback = nullptr ){
    std::shared_ptr<RhsRigidBody> rhs = std::make_shared<RhsRigidBody>(*this);
    std::shared_ptr<NumericDerivative> dlagrange = std::make_shared<NumericDerivative>(rhs);
   
    SolveODE_Alpha (tend, steps, 0.8, q_, dq_, ddq_, dlagrange, mass_function, callback);
  }
};

class RhsRigidBody : public NonlinearFunction
{
  public:
  //Just a reminder
  RhsRigidBody(RigidBody b){};
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
    dNumeric(*this,x,df);
  }

};

