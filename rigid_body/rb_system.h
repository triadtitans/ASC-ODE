#include <vector>
#include "rigid_body.h"
#include <cmath>


class RhsRBSystem;

enum class ConnectorType { mass, fix };

struct Connector{
  ConnectorType t;
  size_t num;
  Vec<3> pos={0,0,0};

  Vec<3> absPos(VectorView<double> x){
    if(t==ConnectorType::fix){
      return pos;
    }
    else{
      Transformation tr (x.Range(num*dim_per_body,(num+1)*dim_per_body));
      Vec<3> res = tr.apply(pos);
      return res;
    }
  }
};

struct Beam{
  double length;
  Connector a;
  Connector b;
};


struct Spring{
  double length;
  double stiffness;
  Connector a;
  Connector b;
};

class RBSystem {
  std::vector<RigidBody> _bodies;
  Vec<3> _gravity={0,0,0};
  std::vector<Beam> _beams;
  std::vector<Spring> _springs;
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
  
  void addSpring(Spring s){
    _springs.push_back(s);
  }

  Connector addFix(){
    return Connector{ConnectorType::fix, 0, {0,0,0}};
  }
  auto& bodies(){return _bodies;}
  auto& gravity(){return _gravity;}
  auto& beams(){return _beams;}
  auto& springs(){return _springs;}
  int numBodies(){return _bodies.size();}
  int numSprings(){return _springs.size();}
  int numBeams(){return _beams.size();}
  int dimension(){return dim_per_body*numBodies()+_beams.size();}

  void GetState(VectorView<double> x, VectorView<double> dx){
    for(int i=0; i<numBodies(); i++){
      x.Range(dim_per_body*i,dim_per_body*i+dim_per_body)=_bodies[i].getQ().q_;
      dx.Range(dim_per_body*i,dim_per_body*i+dim_per_body)=_bodies[i].getDq().q_;
    }
  }
  Vec<3> connectorPos(Connector c){
    Vector<double> x(dimension());
    Vector<double> dx(dimension());
    GetState(x,dx);
    return c.absPos(x);
  }
  
  void SetState(VectorView<double> x, VectorView<double> dx){
    for(int i=0; i<numBodies(); i++){
      _bodies[i].setQ(Transformation(x.Range(dim_per_body*i,dim_per_body*i+dim_per_body)));
      _bodies[i].setDq(Transformation(dx.Range(dim_per_body*i,dim_per_body*i+dim_per_body)));
    }
  }

  void saveState(){
    for(auto& r: bodies())
      r.saveState();
  }

  void reset(){
    for(auto& r: bodies())
      r.reset();
  }
  
  void simulate(double tend, double steps, std::function<void(double,VectorView<double>)> callback = nullptr ){
    //Create Lagragian
    std::shared_ptr<RhsRBSystem> rhs = std::make_shared<RhsRBSystem>(*this);
    //Derive Lagragian
    std::shared_ptr<NumericDerivative> dlagrange = std::make_shared<NumericDerivative>(rhs);
    Vector<double> state(dimension());
    Vector<double> dstate(dimension());

    std::shared_ptr<StackedFunction> mass = std::make_shared<StackedFunction>();
    mass->addFunction(_mass_func);
    //Extend mass func to beam lambdas
    if(_beams.size()>0){
      Vector<double> zero(_beams.size());
      std::shared_ptr<NonlinearFunction> const_zero = std::make_shared<ConstantFunction>(zero);
      mass->addFunction(const_zero);
    }

    GetState(state,dstate);
    SolveODE_Newmark (tend, steps, state, dstate, dlagrange, mass, callback);
    SetState(state,dstate);
    
  } 
};

class RhsRBSystem : public NonlinearFunction
{
  RBSystem& _s;
  public:
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

      //Calculate gravitational potential
      auto t = Transformation(q);
      f(0) -=  _s.bodies()[i].mass()*t.apply(_s.bodies()[i].center())*_s.gravity(); 

    }
    for(int i=0;i<_s.numBeams();i++){
      auto& b = _s.beams()[i];
      //Evaluate distance between connectors, and compare with beam length
      Connector c1 = b.a;
      Connector c2 = b.b;

      Vec<3> pos1 = c1.absPos(x);
      Vec<3> pos2 = c2.absPos(x);
      double g = Norm(pos1-pos2)-b.length;
      f(0)+=g*x(numBodies()*dim_per_body+i);
    }
    for(int i=0;i<_s.numSprings();i++){
      auto& s = _s.springs()[i];
      //Evaluate distance between connectors, and compare with beam length
      Connector c1 = s.a;
      Connector c2 = s.b;

      Vec<3> pos1 = c1.absPos(x);
      Vec<3> pos2 = c2.absPos(x);
      double e = (1/2.0)*s.stiffness*std::pow((Norm(pos1-pos2)-s.length),2);
      f(0)-=e;
    }
  }
  virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
  {
    dNumeric(*this,x,df);
  }

};