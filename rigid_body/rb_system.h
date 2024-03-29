#include <vector>
#include "rigid_body.h"

constexpr int dim_per_body = 18;

class RhsRBSystem;

enum class ConnectorType { mass };

struct Connector{
  ConnectorType t;
  size_t num;
  Vec<3> pos;
};

struct Beam{
  double length;
  Connector a;
  Connector b;
};

class RBSystem {
  std::vector<RigidBody> _bodies;
  std::vector<Beam> _beams;
  std::shared_ptr<StackedFunction> _mass_func;
public:
  RBSystem(){
    _mass_func =std::make_shared<StackedFunction>();
  };
  Connector addBody(RigidBody& b){
    _bodies.push_back(b);
    _mass_func->addFunction(b.getMassFunc());
    return Connector{ConnectorType::mass, _bodies.size()-1, Vector<double>(3)};
  }
  void addBeam(Beam b){
    _beams.push_back(b);
  }
  auto& bodies(){return _bodies;}
  auto& beams(){return _beams;}
  int numBodies(){return _bodies.size();}
  int numBeams(){return _beams.size();}
  int dimension(){return dim_per_body*numBodies()+_beams.size();}

  void GetState(VectorView<double> x, VectorView<double> dx){
    for(int i=0; i<numBodies(); i++){
      x.Range(dim_per_body*i,dim_per_body*i+dim_per_body)=_bodies[i].getQ().q_;
      dx.Range(dim_per_body*i,dim_per_body*i+dim_per_body)=_bodies[i].getDq().q_;
    }
  }

  void SetState(VectorView<double> x, VectorView<double> dx){
    for(int i=0; i<numBodies(); i++){
      _bodies[i].setQ(Transformation(x.Range(dim_per_body*i,dim_per_body*i+dim_per_body)));
      _bodies[i].setDq(Transformation(dx.Range(dim_per_body*i,dim_per_body*i+dim_per_body)));
    }
  }
  
  void simulate(double tend, double steps, std::function<void(double,VectorView<double>)> callback = nullptr ){
    //Create Lagragian
    std::shared_ptr<RhsRBSystem> rhs = std::make_shared<RhsRBSystem>(*this);
    //Derive Lagragian
    std::shared_ptr<NumericDerivative> dlagrange = std::make_shared<NumericDerivative>(rhs);
    Vector<double> state(dimension());
    Vector<double> dstate(dimension());

    //Extend mass func to beam lambdas
    Vector<double> zero(_beams.size());
    std::shared_ptr<NonlinearFunction> const_zero = std::make_shared<ConstantFunction>(zero);
    std::shared_ptr<StackedFunction> mass = std::make_shared<StackedFunction>();
    mass->addFunction(_mass_func);
    mass->addFunction(const_zero);

    GetState(state,dstate);
    SolveODE_Newmark (tend, steps, state, dstate, dlagrange, mass, callback);
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
      Vector<double> q = x.Range(dim_per_body*i,dim_per_body*i+dim_per_body);

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
    for(int i=0;i<_s.numBeams();i++){
      auto& b = _s.beams()[i];
      //Evaluate distance between connectors, and compare with beam length
      Connector c1 = b.a;
      Connector c2 = b.b;

      //Get the position of RigidBody 1 and 2 corresponding to state x
      Transformation t1 (x.Range(c1.num*dim_per_body,(c1.num+1)*dim_per_body));
      Transformation t2 (x.Range(c2.num*dim_per_body,(c2.num+1)*dim_per_body));
      Vec<3> pos1 = t1.apply(c1.pos);
      Vec<3> pos2 = t2.apply(c2.pos);
      double g = Norm(pos1-pos2)-b.length;
      f(0)+=g*x(numBodies()*dim_per_body+i);
    }
  }
  virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
  {
    dNumeric(*this,x,df);
  }

};