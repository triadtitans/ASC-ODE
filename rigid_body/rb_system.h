#include <vector>
#include "rigid_body.h"
#include <cmath>


class RhsRBSystem;

enum class ConnectorType { mass, fix };

struct Connector{
  ConnectorType t;
  size_t num;
  Vec<3> pos={0,0,0};

  template<typename T>
  Vec<3, T> absPos(VectorView<T> x){
    if(t==ConnectorType::fix){
      Vec<3, T> pos_t(pos);
      return pos_t;
    }
    else{
      Transformation<T> tr (x.Range(num*dim_per_body,(num+1)*dim_per_body));
      Vec<3, T> res = tr.apply(pos);
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
      _bodies[i].setQ(Transformation<>(x.Range(dim_per_body*i,dim_per_body*i+dim_per_body)));
      _bodies[i].setDq(Transformation<>(dx.Range(dim_per_body*i,dim_per_body*i+dim_per_body)));
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
    // std::shared_ptr<NumericDerivative> dlagrange = std::make_shared<NumericDerivative>(rhs);
    std::shared_ptr<Derivative> dlagrange = std::make_shared<Derivative>(rhs);
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
  size_t DimX() const  { return _s.dimension();}
  size_t DimF() const  { return 1; }
  int numBodies() const {return _s.numBodies();}

  void Evaluate (VectorView<double> x, VectorView<double> f) const 
  {
    f = 0;
    for(int i=0;i< numBodies(); i++){
      Matrix<AutoDiffDiff<dim_per_body, double>> b (3, 3);

      //x contains state of all bodies
      //q contains state (and lambdas) of one body
      Vector<double> q = x.Range(dim_per_body*i,dim_per_body*i+dim_per_body);

      Vector<AutoDiffDiff<dim_per_body, double>> x_diff(dim_per_body);

      Vector<AutoDiffDiff<dim_per_body, double>> f_diff(1);
    
      for (size_t j = 0; j < dim_per_body; j++) 
      {
        x_diff(j).Value() = q(j);
        x_diff(j).DValue(j) = 1;
      }


   
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

      //Calculate gravitational potential
      auto t = Transformation<AutoDiffDiff<dim_per_body, double>>(x_diff);
      f_diff(0) -=  _s.bodies()[i].mass()*t.apply(_s.bodies()[i].center())*_s.gravity(); 

      for (size_t j = 0; j < dim_per_body; j++) {
        f(dim_per_body*i + j) += f_diff(0).DValue(j);
      }

    }
    for(int i=0;i<_s.numBeams();i++){
      auto& b = _s.beams()[i];
      //Evaluate distance between connectors, and compare with beam length
      Connector c1 = b.a;
      size_t l = c1.num;
      Connector c2 = b.b;
      size_t k = c2.num;

      // set elements of body a
      // set elements of body b
      Vector<AutoDiffDiff<2*dim_per_body + 1, double>> x_diff(DimX());
      for (size_t j = 0; j < dim_per_body; j++) {
        
        x_diff(dim_per_body*l + j).Value() = x(dim_per_body*l + j);
        x_diff(j).DValue(j) = 1;

        x_diff(dim_per_body*k + j).Value() = x(dim_per_body*k + j);
        x_diff(dim_per_body*k + j).DValue(dim_per_body + j) = 1;
      }
      x_diff(numBodies()*dim_per_body+i).Value() = x(numBodies()*dim_per_body+i);
      x_diff(numBodies()*dim_per_body+i).DValue(2*dim_per_body) = 1;

      Vec<3, AutoDiffDiff<2*dim_per_body + 1, double>> pos1 = c1.absPos(x_diff);
      Vec<3, AutoDiffDiff<2*dim_per_body + 1, double>> pos2 = c2.absPos(x_diff);
      //Norm is template (AutoDiffDiff able)
      AutoDiffDiff<2*dim_per_body + 1, double> f_diff = (Norm(pos1-pos2)-b.length)*x_diff(numBodies()*dim_per_body+i);
      
      for (size_t j = 0; j < dim_per_body; j++) {
        f(dim_per_body*l + j) += f_diff.DValue(j);
        f(dim_per_body*k + j) += f_diff.DValue(dim_per_body + j);
      }

      f(numBodies()*dim_per_body+i) += x_diff(numBodies()*dim_per_body+i).DValue(2*dim_per_body);
    }


    for(int i=0;i<_s.numSprings();i++){
      auto& s = _s.springs()[i];
      //Evaluate distance between connectors, and compare with beam length
      Connector c1 = s.a;
      size_t l = c1.num;
      Connector c2 = s.b;
      size_t k = c2.num;

      Vector<AutoDiffDiff<2*dim_per_body, double>> x_diff(DimX());


      for (size_t j = 0; j < dim_per_body; j++) {
        
        x_diff(dim_per_body*l + j).Value() = x(dim_per_body*l + j);
        x_diff(j).DValue(j) = 1;

        x_diff(dim_per_body*k + j).Value() = x(dim_per_body*k + j);
        x_diff(dim_per_body*k + j).DValue(dim_per_body + j) = 1;
      }
      

      Vec<3, AutoDiffDiff<2*dim_per_body, double>> pos1 = c1.absPos(x_diff);
      Vec<3, AutoDiffDiff<2*dim_per_body, double>> pos2 = c2.absPos(x_diff);
      AutoDiffDiff<2*dim_per_body, double> e_diff = (1/2.0)*s.stiffness*pow((Norm(pos1-pos2)-s.length), 2.);
      
      for (size_t j = 0; j < dim_per_body; j++) {
        f(dim_per_body*l + j) += e_diff.DValue(j);
        f(dim_per_body*k + j) += e_diff.DValue(dim_per_body + j);
      }
    }
  }

  
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
  {
    // our function is 2-times 
    df = 0;
    for(int i=0;i< numBodies(); i++){
      Matrix<AutoDiffDiff<dim_per_body, double>> b (3, 3);

      //x contains state of all bodies
      //q contains state (and lambdas) of one body
      Vector<double> q = x.Range(dim_per_body*i,dim_per_body*i+dim_per_body);

      Vector<AutoDiffDiff<dim_per_body, double>> x_diff(dim_per_body);

      Vector<AutoDiffDiff<dim_per_body, double>> f_diff(1);
    
      for (size_t j = 0; j < dim_per_body; j++) 
      {
        x_diff(j).Value() = q(j);
        x_diff(j).DValue(j) = 1;
      }


   
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

      //Calculate gravitational potential
      auto t = Transformation<AutoDiffDiff<dim_per_body, double>>(x_diff);
      f_diff(0) -=  _s.bodies()[i].mass()*t.apply(_s.bodies()[i].center())*_s.gravity(); 

      for (size_t j = 0; j < dim_per_body; j++) {
        for(size_t k = 0; k < dim_per_body; k++) {
          df(dim_per_body*i + j, dim_per_body*i + k) += f_diff(0).DDValue(j, k);
        }
      }

    }
    for(int i=0;i<_s.numBeams();i++){
      auto& b = _s.beams()[i];
      //Evaluate distance between connectors, and compare with beam length
      Connector c1 = b.a;
      size_t l = c1.num;
      Connector c2 = b.b;
      size_t k = c2.num;

      // set elements of body a
      // set elements of body b
      Vector<AutoDiffDiff<2*dim_per_body + 1, double>> x_diff(DimX());
      for (size_t j = 0; j < dim_per_body; j++) {
        
        x_diff(dim_per_body*l + j).Value() = x(dim_per_body*l + j);
        x_diff(j).DValue(j) = 1;

        x_diff(dim_per_body*k + j).Value() = x(dim_per_body*k + j);
        x_diff(dim_per_body*k + j).DValue(dim_per_body + j) = 1;
      }
      x_diff(numBodies()*dim_per_body+i).Value() = x(numBodies()*dim_per_body+i);
      x_diff(numBodies()*dim_per_body+i).DValue(2*dim_per_body) = 1;

      Vec<3, AutoDiffDiff<2*dim_per_body + 1, double>> pos1 = c1.absPos(x_diff);
      Vec<3, AutoDiffDiff<2*dim_per_body + 1, double>> pos2 = c2.absPos(x_diff);
      //Norm is template (AutoDiffDiff able)
      AutoDiffDiff<2*dim_per_body + 1, double> f_diff = (Norm(pos1-pos2)-b.length)*
      x_diff(numBodies()*dim_per_body+i);
      
      for (size_t j = 0; j < dim_per_body; j++) {
        for (size_t h = 0; h < dim_per_body; k++) {
          // diagonalmatritzen
          df(dim_per_body*l + j, dim_per_body*l + h) += f_diff.DDValue(j, h);
          df(dim_per_body*k + j, dim_per_body*k + h) += f_diff.DDValue(dim_per_body + j, dim_per_body + h);
          // Hessematrix elemente body1 gegen body2
          df(dim_per_body*k + j, dim_per_body*l + h) += f_diff.DDValue(dim_per_body + j, h);
          df(dim_per_body*l + j, dim_per_body*k + h) += f_diff.DDValue(j, dim_per_body + h);

        }
        df(numBodies()*dim_per_body+i, dim_per_body*l + j) += f_diff.DDValue(2*dim_per_body, j);
        df(numBodies()*dim_per_body+i, dim_per_body*k + j) += f_diff.DDValue(2*dim_per_body, dim_per_body + j);

        df(dim_per_body*l + j, numBodies()*dim_per_body+i) += f_diff.DDValue(j, 2*dim_per_body);
        df(dim_per_body*k + j, numBodies()*dim_per_body+i) += f_diff.DDValue(dim_per_body + j, 2*dim_per_body);
      }
    }


    for(int i=0;i<_s.numSprings();i++){
      auto& s = _s.springs()[i];
      //Evaluate distance between connectors, and compare with beam length
      Connector c1 = s.a;
      size_t l = c1.num;
      Connector c2 = s.b;
      size_t k = c2.num;

      Vector<AutoDiffDiff<2*dim_per_body, double>> x_diff(DimX());


      for (size_t j = 0; j < dim_per_body; j++) {
        
        x_diff(dim_per_body*l + j).Value() = x(dim_per_body*l + j);
        x_diff(j).DValue(j) = 1;

        x_diff(dim_per_body*k + j).Value() = x(dim_per_body*k + j);
        x_diff(dim_per_body*k + j).DValue(dim_per_body + j) = 1;
      }
      

      Vec<3, AutoDiffDiff<2*dim_per_body, double>> pos1 = c1.absPos(x_diff);
      Vec<3, AutoDiffDiff<2*dim_per_body, double>> pos2 = c2.absPos(x_diff);
      AutoDiffDiff<2*dim_per_body, double> f_diff = (1/2.0)*s.stiffness*pow((Norm(pos1-pos2)-s.length),2.);
      
      for (size_t j = 0; j < dim_per_body; j++) {
        for (size_t h = 0; h < dim_per_body; h++) {
          // diagonalmatritzen
          df(dim_per_body*l + j, dim_per_body*l + h) += f_diff.DDValue(j, h);
          df(dim_per_body*k + j, dim_per_body*k + h) += f_diff.DDValue(dim_per_body + j, dim_per_body + h);
          // Hessematrix elemente body1 gegen body2
          df(dim_per_body*k + j, dim_per_body*l + h) += f_diff.DDValue(dim_per_body + j, h);
          df(dim_per_body*l + j, dim_per_body*k + h) += f_diff.DDValue(j, dim_per_body + h);

        }
      }
    }
  }
};