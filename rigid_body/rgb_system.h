#include <vector>
#include "rigid_body.h"

constexpr int dim_per_body = 18;

class RBSystem {
  std::vector<RigidBody> _bodies;
  std::shared_ptr<StackedFunction> _mass_func;
public:
  RBSystem(){
    _mass_func =std::make_shared<StackedFunction>();
  };
  void addBody(RigidBody& b){
    _bodies.push_back(b);
    _mass_func->addFunction(b.getMassFunc());
  }
  auto& bodies(){return _bodies;}
  int numBodies(){return _bodies.size();}
  int dimension(){return dim_per_body*numBodies();}

  void GetState(VectorView<double> x, VectorView<double> dx){
    for(int i=0; i<numBodies(); i++){
      x.Range(dim_per_body*i,dim_per_body)=_bodies[i].getQ().q_;
      dx.Range(dim_per_body*i,dim_per_body)=_bodies[i].getDq().q_;
    }
  }

  void SetState(VectorView<double> x, VectorView<double> dx){
    for(int i=0; i<numBodies(); i++){
      _bodies[i].setQ(Transformation(x.Range(dim_per_body*i,dim_per_body)));
      _bodies[i].setDq(Transformation(dx.Range(dim_per_body*i,dim_per_body)));
    }
  }
  
  void simulate(double tend, double steps, std::function<void(double,VectorView<double>)> callback = nullptr ){
    std::shared_ptr<RhsRBSystem> rhs = std::make_shared<RhsRBSystem>(*this);
    std::shared_ptr<NumericDerivative> dlagrange = std::make_shared<NumericDerivative>(rhs);
    Vector<double> state(dimension());
    Vector<double> dstate(dimension());
    GetState(state,dstate);
    SolveODE_Newmark (tend, steps, state, dstate, dlagrange, _mass_func, callback);
    SetState(state,dstate);
    
  } 
};

class RhsRBSystem : public NonlinearFunction
{
  RBSystem& _s;
  public:
  //Just a reminder
  RhsRBSystem(RBSystem& s):_s(s){};
  size_t DimX() const override { return _s.dimension(); }
  size_t DimF() const override { return 1; }
  int numBodies() const {return _s.numBodies();}
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0)=0;
    for(int i=0;i< numBodies(); i++){
      Matrix<double> b (3, 3);

      //x contains state of all bodies
      //q contains state (and lambdas) of one body
      Vector<double> q = x.Range(dim_per_body*i,dim_per_body);

      // extract B from Q
      b.Row(0) = q.Range(1, 4);
      b.Row(1) = q.Range(5, 8);
      b.Row(2) = q.Range(9, 12);
      // b must be orthonormal
      auto c = TransposeMatExpr(b)* b + (-1)*IdMatExpr(3);
      Vector<double> g(6);
      g(0)=c(0, 0);
      g(1)=c(1, 0);
      g(2)=c(2, 0);
      g(3)=c(1, 1);
      g(4)=c(2, 1);
      g(5)=c(2, 2);
      f(0) += q.Range(12, 18)*g;
    }
  }
  virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
  {
    dNumeric(*this,x,df);
  }

};