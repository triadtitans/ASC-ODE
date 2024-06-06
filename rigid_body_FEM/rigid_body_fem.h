#include <nonlinfunc.h>
#include <math.h>
#include <matrix.h>
#include <vector.h>
#include <ode.h>
#include "../src/autodiffdiff.h"
#include <cmath>


#ifdef PYBIND11_MODULE
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif


constexpr size_t dim_per_transform = 12;
constexpr size_t dim_per_body = 18;
constexpr size_t eq_per_body = 30;

class RigidBody_FEM;
class RBS_FEM;
class EQRigidBody;
class EQRBSystem;

using namespace ASC_bla;
using namespace ASC_ode;


// convert vector v to skew-symmetric matrix
template<typename T>
auto hat_map(const VecExpr<T>& v){
  Matrix<decltype(1.0*v(0))> M(3, 3);
  M(0, 1) = -v(2);
  M(0, 2) = v(1);
  M(1, 0) = v(2);
  M(1, 2) = -v(0);
  M(2, 0) = -v(1);
  M(2, 1) = v(0);
  return M;
}

// convert skew-symmetric matrix to vector
template<typename T>
auto inv_hat_map(const MatrixExpr<T>& M){
  Vector<decltype(1.0 * M(0, 0))> v(3);
  v(0) = (M(2, 1) - M(1, 2))/2;
  v(1) = (M(0, 2) - M(2, 0))/2;
  v(2) = (M(1, 0) - M(0, 1))/2;
  return v;
}


template<typename T = double>
class Transformation{
 public:
  // ordering of q_: Schöberl-style
  Vector<T> q_;

  Transformation(Vector<T> q):q_(q){}
  Transformation():q_(12){
    q_(1)=1;
    q_(6)=1;
    q_(11)=1;
  }

  Vec<3, T> apply(Vec<3, T> pos){
    Vec<3, T> trans{q_(0),q_(4),q_(8)};
    Matrix<T> rot = AsMatrix(q_,3,4).Cols(1,3);
    return trans + rot*pos;
  }

  void setTranslation(T a, T b, T c){q_(0)=a;q_(4)=b;q_(8)=c;}
  void setRotation(int i, int j, T r){
    if(i>2||j>2 || j<0 || i <0) throw std::invalid_argument("Rotation Matrix is 3x3");
    q_((4*i + 1) + j)=r; // note the ordering of q
  }
  template<typename S>
  void setRotation_from_matrix(const MatrixExpr<S>& B){
    for (size_t i=0; i < 3; i++){
      for (size_t j=0; j < 3; j++){
        setRotation(i, j, B(i, j));
      }
    }
  }
  // set rotation matrix from x, y or z degrees
  void setRotationDeg(int axis, T deg){
    Matrix<T> R = makeRotationMatrix3<T>(axis, deg);
    setRotation_from_matrix(R);
  }
  // add a rotation around axis (by index) with degrees deg
  void applyRotationDeg(int axis, T deg){
    Matrix<T> R_new = makeRotationMatrix3<T>(axis, deg);
    Matrix<T> R_old = getRotation();
    setRotation_from_matrix(R_old*R_new);
  }

  Vector<T> getTranslation() const{
    Vector<T> v(3);
    v(0)=q_(0);
    v(1)=q_(4);
    v(2)=q_(8);
    return v;
  }
  Matrix<T> getRotation() const {
    Matrix<T> B(3, 3);
    for (size_t i=0; i < 3; i++){
      for (size_t j=0; j < 3; j++){
        B(i, j) = q_((4*i + 1) + j);
      }
    }
    return B;
  }
};

  std::ostream& operator<<(std::ostream& oss, const Transformation<double>& t){
      oss<<std::fixed << " Translation: \t" << t.q_(0) << " ," << t.q_(4) << ", "<<", " << t.q_(8) << std::endl
                      << " Rotation: \t" << t.q_(1) << " ," << t.q_(2) << ", "<<", " << t.q_(3) <<  std::endl
                     << " \t\t" << t.q_(5) << " ," << t.q_(6) << ", "<<", " << t.q_(7) <<  std::endl
                      << "\t\t" << t.q_(9) << " ," << t.q_(10) << ", "<<", " << t.q_(11) << std::endl; 
    return oss;
  }



enum class ConnectorType { mass, fix };
// manages connections to bodies with Springs etc.
struct Connector{
  ConnectorType t;
  // size_t num;
  Vec<3> pos={0,0,0};
  size_t body_index; // index of body in system's _bodies that the Connector is relative to; irrelevant for fixes

  template<typename T>
  Vec<3, T> absPos(Vector<T> a, Matrix<T> B) const{
    if(t == ConnectorType::fix){
      Vec<3, T> pos_t(pos);
      return pos_t;
    }
    else{
      Transformation<T> tr;
      tr.setTranslation(a(0), a(1), a(2));
      tr.setRotation_from_matrix(B);
      Vec<3, T> res = tr.apply(pos);
      return res;
    }
  }
};

struct Spring{
  double length;
  double stiffness;
  Connector a;
  Connector b;
};

struct Beam{
  double length;
  Connector a;
  Connector b;
};


class RigidBody_FEM {
  Vector<double> q_;
  Vector<double> phat_;
  Vector<double> initialq_;
  Vector<double> initialphat_;

  Matrix<double> P_;

  Vec<3> center_of_mass_={0,0,0};
  Matrix<double> inertia_;
  double mass_=1;
  std::shared_ptr<NonlinearFunction> mass_function;
  Matrix<double> mass_matrix_;

  #ifdef PYBIND11_MODULE
  py::list vertices_;
  py::list normals_;
  #endif

public:

  template<typename T>
  RigidBody_FEM(Vector<double> q,Vector<double> phat,double mass, Vec<3> center_of_mass, MatrixExpr<T>& inertia)
        : q_(q), phat_(phat), initialq_(q), initialphat_(phat), center_of_mass_(center_of_mass),
          mass_(mass), inertia_(inertia), P_(6, 12), mass_matrix_(Diagonal<double>(6, 1.0)) {
    if(inertia.Width() != 3) throw std::invalid_argument("Inertia matrix must be 3x3");
    if(inertia.Height() != 3) throw std::invalid_argument("Inertia matrix must be 3x3");
    //if(q.Size() != dim_per_body) throw std::invalid_argument("q Vector must match mass matrix");
    if(phat.Size() != 6) throw std::invalid_argument("q Vector must match mass matrix");
  }

  RigidBody_FEM()
        :   mass_matrix_(Diagonal<double>(6, 1.0)),
          q_(12), initialq_(12), phat_(6), initialphat_(6), inertia_(3,3),center_of_mass_{0,0,0}, P_(6, 12){
    q_(1)=1;q_(6)=1;q_(11)=1;
    // set inertia?
  };
  
 
  double& mass(){return mass_;}
  Vec<3>& center(){return center_of_mass_;}
  Matrix<double>& inertia(){return inertia_;}

  #ifdef PYBIND11_MODULE
  py::list& vertices(){return vertices_;}
  py::list& normals(){return normals_;}

  py::tuple to_pickle(){ // __getstate__
    auto vectorToTuple = [](Vector<double>& v){
      return py::make_tuple(v.Size(),
                            py::bytes((char*)(void*)&v(0), v.Size()*sizeof(double)));
    };
    auto matrixToTuple = [](Matrix<double>& m){
      return py::make_tuple(m.Height(), m.Width(),
                            py::bytes((char*)(void*)&m(0,0), m.Height()*m.Width()*sizeof(double)));
    };
    py::tuple t_q = vectorToTuple(q_);
    py::tuple t_phat = vectorToTuple(phat_);
    py::tuple t_center_of_mass = py::make_tuple(center_of_mass_(0),center_of_mass_(1),center_of_mass_(2));
    py::tuple t_inertia = matrixToTuple(inertia_);

    return py::make_tuple(t_q, t_phat, t_center_of_mass, t_inertia, mass_, vertices_, normals_);
  }
  void load_pickle(py::tuple data){ // __setstate__
    if (data.size() != 7) throw invalid_argument("invalid tuple length for pickled object");

    auto tupleToVector = [](py::tuple t){ // __setstate__
      if (t.size() != 2)
        throw std::runtime_error("should be a 2-tuple!");

      Vector<double> v(t[0].cast<size_t>());
      py::bytes mem = t[1].cast<py::bytes>();
      std::memcpy(&v(0), PYBIND11_BYTES_AS_STRING(mem.ptr()), v.Size()*sizeof(double));
      return v;
    };
    auto tupleToMatrix = [](py::tuple t){ // __setstate__
      if (t.size() != 3)
        throw std::runtime_error("should be a 3-tuple!");

      Matrix<double> m(t[0].cast<size_t>(),t[1].cast<size_t>());
      py::bytes mem = t[2].cast<py::bytes>();
      std::memcpy(&m(0,0), PYBIND11_BYTES_AS_STRING(mem.ptr()), m.Height()*m.Width()*sizeof(double));
      return m;
    };

    q_ = tupleToVector(data[0]);
    phat_ = tupleToVector(data[1]);
    initialq_ = tupleToVector(data[0]);
    initialphat_ = tupleToVector(data[1]);
    py::tuple ecom = data[2]; // entry of center of mass
    center_of_mass_(0) = (ecom[0]).cast<double>();
    center_of_mass_(1) = (ecom[1]).cast<double>();
    center_of_mass_(2) = (ecom[2]).cast<double>();
    inertia_ = tupleToMatrix(data[3]);
    mass_ = data[4].cast<double>();
    recalcMassMatrix();
    vertices_ = data[5];
    normals_ = data[6];
  }
  #endif

  void recalcMassMatrix(){
    mass_matrix_ = Matrix(6, 6);
    mass_matrix_(0, 0) = 1.0;
    mass_matrix_(1, 1) = 1.0;
    mass_matrix_(2, 2) = 1.0;
    mass_matrix_ = mass_*mass_matrix_;
    mass_matrix_.Rows(3, 3).Cols(3, 3) = inertia_; // inertia matrix is already multiplied with mass
  }

  Vector<double>& q(){return q_;}
  Vector<double>& phat(){return phat_;}
  Matrix<double>& Mass_matrix(){return mass_matrix_;}

  void setQ(Transformation<> t){q_=t.q_;}
  void setPhat(Vector<double> v){phat_=v;}
  void setPhat_v(size_t i, double a) {
    phat_(i) = a;
  }

  void setMass(Matrix<double> m){mass_function=std::make_shared<LinearFunction>(m);}
  std::shared_ptr<NonlinearFunction> getMassFunc(){return mass_function;}
  // saves a state for the reset button
  void saveState(){
    initialq_ = q_;
    initialphat_ = phat_;
  }
  // resets position and rotation to last saved state
  void reset(){
    q_ = initialq_;
    phat_ = initialphat_;
  }

  Vec<3> absolutePosOf(Vec<3> relative_pos){
    // std::cout << getQ() << std::endl;
    return getQ().apply(relative_pos);
  }

  Transformation<> getQ(){return q_;}
  Vector<double> getPhat(){return phat_;}

};



template<typename T>
Vector<T> get_translation(VectorView<T> x, size_t body_index)  {
  if(x.Size()%dim_per_transform)
    throw std::invalid_argument("Vector must be in transform format");

  Vector<T> trafo = {x(body_index*dim_per_transform + 0), x(body_index*dim_per_transform + 1), x(body_index*dim_per_transform + 2)};
  return trafo;
}

template<typename T>
Matrix<T> get_rotation(VectorView<T> x, size_t body_index)  {
  if(x.Size()%dim_per_transform)
    throw std::invalid_argument("Vector must be in transform format");
  
  Matrix<T> B = AsMatrix(x.Range(body_index*dim_per_transform + 3, body_index*dim_per_transform + 12), 3, 3);

  return B;
}

class RBS_FEM{
  std::vector<RigidBody_FEM> _bodies;
  std::vector<Spring> _springs;
  std::vector<Beam> _beams;
  Vector<double> gravity_ = {0, 0, 0};

  public:
  std::vector<RigidBody_FEM>& bodies(){return _bodies;}
  std::vector<Spring>& springs(){return _springs;}
  auto& beams(){return _beams;}
  Vector<double> & gravity() {return gravity_;}
  int numBodies(){return _bodies.size();}
  int numSprings(){return _springs.size();}
  int numBeams(){return _beams.size();}

  void getState(VectorView<double> out){
    for(int i=0; i<_bodies.size(); i++){
      out.Range(i*dim_per_body,i*dim_per_body+12)=_bodies[i].q();
      out.Range(i*dim_per_body+12, (i+1)*dim_per_body)=_bodies[i].phat();
    }
  }
  void setState(VectorView<double> in){
    for(int i=0; i<_bodies.size(); i++){
      _bodies[i].q()=in.Range(i*dim_per_body,i*dim_per_body+12);
      _bodies[i].phat()=in.Range(i*dim_per_body+12, (i+1)*dim_per_body);
    }
  }

  // convert (Newton's) equation solution format to rigid body system state format
  // output format: q, phat, q, phat, q, phat, ... for the different bodies in _bodies
  // ordering of a single q in output: as stored by Transformation
  // input format: see EQRigidBody::Evaluate
  Vector<double> xToState(VectorView<double> x){
    if(x.Size()%eq_per_body)
      throw std::invalid_argument("Vector must be in Equation format");

    size_t num_bodies=x.Size()/eq_per_body;
    // the Vector in output format:
    Vector<double> res(num_bodies * dim_per_body);

    for(size_t i = 0;i<num_bodies; i++){
      // make Transformation object from equation solution format
      Transformation<double> T;
      T.setTranslation(x(i*eq_per_body + 0), x(i*eq_per_body + 1), x(i*eq_per_body + 2));
      T.setRotation_from_matrix(AsMatrix(x.Range(i*eq_per_body + 3, i*eq_per_body + 12), 3, 3));
      
      // store transformation data into array
      res.Range(i*dim_per_body,i*dim_per_body+12)=T.q_;
      res.Range(i*dim_per_body+12,i*dim_per_body+18)=x.Range(i*eq_per_body+18,i*eq_per_body+24);
    }
    return res;
  }

  Vector<double> stateToX (VectorView<double>x ){
    if(x.Size()%dim_per_body)
      throw std::invalid_argument("Vector must be in Equation format");
    size_t num_bodies=x.Size()/dim_per_body ;
    Vector<double> res(num_bodies * eq_per_body);
    for(size_t i = 0;i<num_bodies; i++){
      Transformation<double> Q;
      Q.q_=x.Range(i*dim_per_body,i*dim_per_body+12);

      res.Range(i*eq_per_body,i*eq_per_body+3)=Q.getTranslation();
      AsMatrix(res.Range(i*eq_per_body+3,i*eq_per_body+12), 3, 3) = Q.getRotation();
      

      res.Range(i*eq_per_body+18,i*eq_per_body+24)=x.Range(i*dim_per_body+12,i*dim_per_body+18);
    }
    return res;
  }

   Vector<double> xToTransformations(VectorView<double> x){
    if(x.Size()%eq_per_body)
      throw std::invalid_argument("Vector must be in Equation format");

    size_t num_bodies=x.Size()/eq_per_body;
    // the Vector in output format:
    Vector<double> res(num_bodies * dim_per_transform);

    for(size_t i = 0;i<num_bodies; i++){
      // store transformation data into array
      res.Range(i*dim_per_transform,i*dim_per_transform+dim_per_transform)=x.Range(i*eq_per_body, i*eq_per_body + dim_per_transform);
      // res.Range(i*dim_per_body+12,i*dim_per_body+18)=x.Range(i*eq_per_body+18,i*eq_per_body+24);
    }
    return res;
  }

  void saveState(){
    for (auto& rb: _bodies){
      rb.saveState();
    }
  }

  void reset(){
    for(auto& body: bodies())
      body.reset();
  }

  Connector addBody(RigidBody_FEM& b){
    _bodies.push_back(b);
    // std::cout << _bodies.size()-1 << std::endl;
    return Connector{ConnectorType::mass, Vector<double>(3), _bodies.size()-1};
  }
  Vec<3> connectorPos(Connector c){
    /* Vector<double> state(18);
    getState(state); */
    /* Vector<double> a = ;
    Matrix<double> B = ;
    return c.absPos(a, B); */
    Vec<3> pos;
    if (c.t == ConnectorType::mass){
      // c.pos is relative position
      return _bodies[c.body_index].absolutePosOf(c.pos);
    }
    else if (c.t == ConnectorType::fix){
      // c.pos is already absolute
      return c.pos;
    }
    else throw invalid_argument("unknown ConnectorType");
    // std::cout << "cp" << relpos << std::endl << _bodies[c.body_index].absolutePosOf(relpos) << std::endl;
    // std::cout << c.body_index << std::endl;
  }

  void addSpring(Spring s){
    _springs.push_back(s);
  }

  void addBeam(Beam b){
    _beams.push_back(b);
  }

  Connector addFix(){
    return Connector{ConnectorType::fix, {0,0,0}, 0};
  }

};



// system of equations for one rigid body
class EQRigidBody : public NonlinearFunction
{
  double h_; // timestep duration
  Transformation<double> Q_; // current body transformation
  Vector<double> phatold; // current momentum
  RBS_FEM& rbs_; // parent rigid body system
  size_t body_index_; // index within _bodies of parent rbs_

  Vector<double> force_half_;
  Vector<double> force_new_;
  
  public:
  EQRigidBody(Transformation<double> Q, Vector<double> Phat, double h, RBS_FEM& rbs, size_t body_index):
      h_(h), Q_(Q), phatold(Phat), rbs_(rbs), body_index_(body_index), force_half_(dim_per_transform), force_new_(dim_per_transform) {};
  size_t DimX() const  { return 30; }
  size_t DimF() const  { return 30; }

  Vector<double>& force_half(){
    return force_half_;
  }
  Vector<double>& force_new(){
    return force_new_;
  }

  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    // x (a(3), B(row_maj)(9), v_trans(3), v_skew(3), phat(6), p(6))

    // variables
    Vector<double> anew = x.Range(0, 3);
    MatrixView<double> Bnew = AsMatrix(x.Range(3, 12), 3, 3);
    Vector<double> vtrans = x.Range(12, 15);
    Vector<double> vskew = x.Range(15, 18);
    Vector<double> phat = x.Range(18, 24);
    Vector<double> p = x.Range(24, 30);

    // known constants
    Vector<double> aold = Q_.getTranslation();
    Matrix<double> Bold = Q_.getRotation();

    Matrix<double> Bhalf(3, 3);
    Bhalf = 0.5*(Bnew + Bold);

    // std::cout << vskew << "\n" << std::endl;

    // I
    f.Range(0, 3) = (1/h_)*(anew - aold) - vtrans;
    f.Range(3, 6) = (1/h_)*inv_hat_map(Transpose(Bhalf)*(Bnew - Bold)) - vskew;

    // II
    Vector<double> vhat(6);
    vhat.Range(0, 3) = vtrans;
    vhat.Range(3, 6) = vskew;

    // Mass matrix included
    f.Range(6, 12) = rbs_.bodies()[body_index_].Mass_matrix()*vhat - p;

    //TODO: Are Bold and Bnew in the right order?
    // III - first half
    /* Vector<double> a_force(3);
    Matrix<double> B_force(3, 3); */
    
    // force<double, double> (0.5*(aold+anew), Bhalf, a_force, B_force); // half or old?
    f.Range(12, 15) = (2/h_)*(p.Range(0, 3) - phatold.Range(0, 3)) - get_translation(force_half_, 0);
    f.Range(15, 18) = inv_hat_map((Transpose(Bold) * ((2/h_)*(Bhalf*hat_map(p.Range(3, 6)) - Bold*hat_map(phatold.Range(3, 6))) - get_rotation(force_half_, 0)))); //TODO: Are Bold and Bnew in the right order?

    // III - second half
    // force<double, double> (anew, Bnew, a_force, B_force);
    f.Range(18, 21) = (2/h_)*(phat.Range(0, 3) - p.Range(0, 3)) - get_translation(force_new_, 0);
    // std::cout << get_translation(force_new_, 0) << std::endl;
    f.Range(21, 24) = inv_hat_map((Transpose(Bnew) * ((2/h_)*(Bnew*hat_map(phat.Range(3, 6)) - Bhalf*hat_map(p.Range(3, 6))) - get_rotation(force_new_, 0)))); //TODO: Are Bold and Bnew in the right order?

    // Bnew must be orthonormal
    Matrix<double> eye = Diagonal(3, 1);
    auto c = Transpose(Bnew)* Bnew + (-1)*eye;
    f(0 + 24)=c(0, 0);
    f(1 + 24)=c(1, 0);
    f(2 + 24)=c(2, 0);
    f(3 + 24)=c(1, 1);
    f(4 + 24)=c(2, 1);
    f(5 + 24)=c(2, 2);
  }

  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    dNumeric(*this, x, df);
  }
};

class EQRigidBodySystem : public NonlinearFunction
{
  RBS_FEM& rbs_;
  std::shared_ptr<StackedFunction> _func;
  std::vector<std::shared_ptr<EQRigidBody>> _functions;
  size_t _num_bodies;

 public:
  size_t DimX() const  { return eq_per_body*_num_bodies; }
  size_t DimF() const  { return  eq_per_body*_num_bodies; }

  EQRigidBodySystem(RBS_FEM& rbs, Vector<double> state, size_t num_bodies, double h): rbs_(rbs), _num_bodies(num_bodies){
    _func = std::make_shared<StackedFunction>();
    // rechnen potential und force => welche parameter für 1 body =>
    for(int i=0;i<num_bodies;i++){
      auto q = state.Range(i*dim_per_body,i*dim_per_body+12);
      auto p = state.Range(i*dim_per_body+12,i*dim_per_body+18); 
      std::shared_ptr<EQRigidBody> eq = std::make_shared<EQRigidBody>(Transformation<double>(q), p, h, rbs_, i); //(Transformation<double>(q), p, h, rbs_, i);
      _func->addFunction(eq);
      _functions.push_back(eq);
    }
  }




  template <typename T>
  T calculate_gravitation_V(Vector<T> a, Matrix<T> B, size_t body_index) const {
    auto t = Transformation<T>();
    t.setRotation_from_matrix(B);
    t.setTranslation(a(0), a(1), a(2));
    
    return (-1)*rbs_.bodies()[body_index].mass() * t.apply(rbs_.bodies()[body_index].center()) * rbs_.gravity();
  }

  
  template<typename T>
  T V (VectorView<T> transforms) const {
    T potential = 0;

    // gravity
    for (size_t i=0; i < _num_bodies; i++){
      potential += calculate_gravitation_V(get_translation(transforms, i), get_rotation(transforms, i), i);
    }

    // spring potential
    for (Spring& spr: rbs_.springs()) {
      // auto& s = rbs_.springs()[i];
      //Evaluate distance between connectors, and compare with beam length
      Connector c1 = spr.a;
      size_t l = c1.body_index;
      Connector c2 = spr.b;
      size_t k = c2.body_index;

      Vec<3, T> pos1 = c1.absPos(get_translation(transforms, l), get_rotation(transforms, l));
      Vec<3, T> pos2 = c2.absPos(get_translation(transforms, k), get_rotation(transforms, k));
      T norm = Norm(pos1-pos2)-spr.length;
      potential += (1/2.0)*spr.stiffness*(norm * norm);
    }

    return potential;
  }

  template<typename T, typename S>
  void force(VectorView<T> transform, VectorView<S>& f) const {
    Vector<AutoDiffDiff<1, T>> transform_diff = transform;
    
    for (size_t i=0; i < transform.Size(); i++){
      transform_diff(i).DValue(0) = 1;
      f(i) = V<>(transform_diff).DValue(0);
      // std::cout << V<>(transform_diff).Value() << std::endl;
      transform_diff(i).DValue(0) = 0;
    }

    /* double eps = 1e-3;
    Vector<> xl(transform.Size()), xr(transform.Size()), fl(1), fr(1);
    for (size_t i = 0; i < transform.Size(); i++)
      {
        xl = transform;
        xl(i) -= eps;
        xr = transform;
        xr(i) += eps;
        fl = V (xl);
        fr = V (xr);
        f(i) = 1/(2*eps) * (fr(0)-fl(0));
      } */
  }

  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    
    Vector<double> transformations_new = rbs_.xToTransformations(x);

    // calculate forces globally
    Vector<double> forces_half(transformations_new.Size());
    Vector<double> state_old(dim_per_body*_num_bodies);
    rbs_.getState(state_old);
    Vector<double> transformations_half = 0.5*(rbs_.xToTransformations(rbs_.stateToX(state_old)) + transformations_new);
    force<>(transformations_half, forces_half);
    Vector<double> forces_new(transformations_new.Size());
    force<>(transformations_new, forces_new);

    // std::cout << "new: " << forces_new << std::endl << "half: " << forces_half << std::endl;
    // std::cout << "new: " << transformations_new << std::endl << "half: " << transformations_half << std::endl;
    // std::cout << _num_bodies << std::endl;
    for (size_t i = 0; i < _num_bodies; i++){
      _functions[i]->force_half() = forces_half.Range(i * dim_per_transform, (i+1) * dim_per_transform);
      _functions[i]->force_new() = forces_new.Range(i * dim_per_transform, (i+1) * dim_per_transform);
      // std::cout << "new: " << _functions[i]->force_new() << std::endl << "half: " << _functions[i]->force_half() << std::endl;
    }
    _func->Evaluate(x,f);
  }
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    _func->EvaluateDeriv(x,df);
  }
};

/* void simulateRigidBody(RigidBody_FEM& rb, double tend, double steps, std::function<void(int,double,VectorView<double>)> callback = nullptr ){
  // solution variable for newton
  Vector<double> state(30);
  // copy over current values 
  state.Range(0,3)=rb.getQ().getTranslation();
  AsMatrix(state.Range(3, 12), 3, 3) = rb.getQ().getRotation();
  state.Range(18, 24) = rb.getPhat();

  for (size_t step=0; step < steps; step++){
    // set up equation
    auto q = rb.getQ();
    auto p = rb.getPhat();
    std::shared_ptr<EQRigidBody> eq = std::make_shared<EQRigidBody>(q, p, tend/steps);
    Matrix<double> df(30,30);
    eq->EvaluateDeriv(state,df);
    std::cout << df.Rows(0,30).Cols(0,30) << std::endl;
    Vector<double> test(30);
      eq->Evaluate(state,test);
          std::cout << test<< std::endl;
    // solve equation
    NewtonSolver(eq, state, 1e-10, 10, callback);

    //store data
    Transformation<double> Q = rb.getQ();
    Q.setTranslation(state(0),state(1),state(2));
    Q.setRotation_from_matrix(AsMatrix(state.Range(3, 12), 3, 3));
    rb.setQ(Q);
    rb.phat() = state.Range(18, 24);
  }
} */

void simulate(RBS_FEM& rbs, double tend, double steps, std::function<void(int,double,VectorView<double>)> callback = nullptr ){  

  Vector<double> state(dim_per_body*rbs.bodies().size());
  Vector<double> x(eq_per_body*rbs.bodies().size());

  // copy over current values 
  rbs.getState(state);
  // std::cout << state << std::endl;
  x = rbs.stateToX(state);
  // std::cout << x<< std::endl;
  for (size_t step=0; step < steps; step++){
    std::shared_ptr<EQRigidBodySystem> eq = std::make_shared<EQRigidBodySystem>(rbs, state,rbs.bodies().size(), tend/steps);
    // solve equation
    Matrix<double> df(60,60);
    //Vector<double> test(30);
    //eq->Evaluate(x,test);
    //auto mv = df.Cols(0,30).Rows(0,30);
    //Matrix<double> mv(30,30);
    //eq->EvaluateDeriv(x,df);
    //std::cout << test<< std::endl;
    //std::cout << df << std::endl;
    // std::cout << "before newton" << std::endl;
    NewtonSolver(eq, x, 1e-10, 10, callback);

    state = rbs.xToState(x);

    // std::cout << x << "\n" << std::endl;
    // std::cout << state << std::endl;

    //store data into different bodies
    rbs.setState(state);

    /* Transformation<double> t = rbs.bodies()[1].q();
    std::cout << t << std::endl; */

  }

} 
