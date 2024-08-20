#ifndef rbs_clean
#define rbs_clean
#include <nonlinfunc.h>
#include <math.h>
#include <matrix.h>
#include <vector.h>
#include <ode.h>
#include "../src/autodiffdiff.h"
#include "rigid_body_helper.h"
#include <cmath>
#include <chrono>


class RgigBody_FEM;
class RBS_FEM;

using namespace ASC_bla;
using namespace ASC_ode;


class RigidBody_FEM {
  //  postion at the end of a timestep
  Vector<double> q_;
  //  momentum at end of a timestep
  Vector<double> phat_;
  //  force at the end of a timestep
  Vector<double> force_;
  //  derivative of the velocity constraint at the end of a time step
  Vector<double> vel_con_;
  //  position at the beginning of the simulation
  Vector<double> initialq_;
  //  momentum at the beginning of the simulation
  Vector<double> initialphat_;
  //  momentum at the beginning of the 

  Matrix<double> P_;

  //  the center of the mass in coordinates relative to body psition
  Vec<3> center_of_mass_={0,0,0};

  //  interia matrix of body (is being set via python and occ)
  Matrix<double> inertia_;

  //  mass of the object if it has no changing mass in the body
  double mass_=1;

  //  mass function - gives the mass at every point in the body
  std::shared_ptr<NonlinearFunction> mass_function_;

  //  mass matrix of body
  Matrix<double> mass_matrix_;

  Matrix<double> mass_matrix_large_;

  //  the index of the body in the system; -1 if the body is not added to a system
  size_t index_ = -1;

  //  list of springs connected to body
  std::vector<size_t> springs_{};

  //  list of beams connected to body
  std::vector<size_t> beams_;

public:

  template<typename T>
  RigidBody_FEM(Vector<double> q,Vector<double> phat,double mass, Vec<3> center_of_mass, MatrixExpr<T>& inertia)
        : q_(q), phat_(phat), force_(Vector<double>(dim_per_transform)), vel_con_(Vector<double>(dim_per_state)),
         initialq_(q), initialphat_(phat), center_of_mass_(center_of_mass), mass_(mass), inertia_(inertia), P_(6, 12), 
         mass_matrix_(Diagonal<double>(6, 1.0)), mass_matrix_large_(Diagonal<double>(3, 1.0)) {
    if(inertia.Width() != 3) throw std::invalid_argument("Inertia matrix must be 3x3");
    if(inertia.Height() != 3) throw std::invalid_argument("Inertia matrix must be 3x3");
    //if(q.Size() != dim_per_body) throw std::invalid_argument("q Vector must match mass matrix");
    if(phat.Size() != 6) throw std::invalid_argument("q Vector must match mass matrix");
  }

  RigidBody_FEM()
        :   mass_matrix_(Diagonal<double>(6, 1.0)),
          q_(12), initialq_(12), phat_(6), initialphat_(6), inertia_(3,3),center_of_mass_{0,0,0}, P_(6, 12), 
          force_(Vector<double>(dim_per_transform)), vel_con_(Vector<double>(dim_per_state)), mass_matrix_large_(Diagonal<double>(3, 1.0)) {
    q_(1)=1;q_(6)=1;q_(11)=1;
    // set inertia?
  };

  size_t& Index()  {
    return index_;
  }
  std::vector<size_t>& Springs()  {
    return springs_;
  }
  std::vector<size_t> Springs() const {
    return springs_;
  }
  size_t NumSprings() const {
    return springs_.size();
  }
  std::vector<size_t>& Beams() {
    return beams_;
  }
  std::vector<size_t> Beams() const {
    return beams_;
  }
  size_t NumBeams() const {
    return beams_.size();
  }
  double& mass()  {
    return mass_;
  }
  double mass() const {
    return mass_;
  }
  Vec<3>& center()  {
    return center_of_mass_;
  }
  Vec<3> center()  const  {
    return center_of_mass_;
  }
  Vector<double>& Force()  {
    return force_;
  }
  Vector<double>& VelocityConstraint()  {
    return vel_con_;
  }
  Matrix<double>& inertia() {
    return inertia_;
  }
  Vector<double>& q() {
    return q_;
  }
  Vector<double>& phat()  {
    return phat_;
  }
  Matrix<double>& Mass_matrix() {
    return mass_matrix_;
  }
  Matrix<double>& Mass_matrix_large_invers()  {
    return mass_matrix_large_;
  }
  void setQ(Transformation<> t) {
    q_=t.q_;
  }
  void setPhat(Vector<double> v)  {
    phat_=v;
  }
  void setPhat_v(size_t i, double a) {
    phat_(i) = a;
  }

  void recalcMassMatrix_large_inverse(){
    //mass_matrix_large_(0,0) = mass_;
    //mass_matrix_large_(1,1) = mass_;
    //mass_matrix_large_(2,2) = mass_;

    Matrix<double> diagblock(3,3);
    for (size_t i=1; i < 4; i++)  {
    diagblock(i, 0) = mass_*center()(i-1);
    diagblock(0, i) = mass_*center()(i-1);
    }
  
    // notes in Norbert's Theory/Massenmatrix.xopp
    // lower right corner (off-diagonal) :
    for (size_t i=1; i < 4; i++){
      for (size_t j=1; j < 4; j++){
        if (i != j){
          diagblock(i, j) = inertia_(i-1, j-1);
        }
      }
    }

    // the diagonal of the lower right corner:
    diagblock(0+1, 0+1) = (inertia_(1, 1) + inertia_(2, 2) - inertia_(0, 0))/2.;
    diagblock(1+1, 1+1) = (inertia_(0, 0) + inertia_(2, 2) - inertia_(1, 1))/2.;
    diagblock(2+1, 2+1) = (inertia_(0, 0) + inertia_(1, 1) - inertia_(2, 2))/2.;

    //mass_matrix_large_.Rows(3,6).Cols(3,6) = diagblock;
    //mass_matrix_large_.Rows(6,9).Cols(6,9) = diagblock;
    //mass_matrix_large_.Rows(9,12).Cols(9,12) = diagblock;
    std::cout << diagblock << std::endl;

    mass_matrix_large_ = inverse(diagblock);
  }

  //  calculates a specific format of the mass matrix fitting to the variable notation of our simulation
  //  Inertia (3x3) 0-Block (3x3)
  //  0-Block (3x3) 
  void recalcMassMatrix() {
    mass_matrix_ = Matrix(6, 6);
    mass_matrix_(0, 0) = 1.0;
    mass_matrix_(1, 1) = 1.0;
    mass_matrix_(2, 2) = 1.0;
    mass_matrix_ = mass_*mass_matrix_;
    mass_matrix_.Rows(3, 3).Cols(3, 3) = inertia_; // inertia matrix is already multiplied with mass
  }

  void setMass(Matrix<double> m)  {
    mass_function_=std::make_shared<LinearFunction>(m);
  }

  std::shared_ptr<NonlinearFunction> getMassFunc()  {
    return mass_function_;
  }
  
  // saves a state for the reset button
  void saveState()  {
    initialq_ = q_;
    initialphat_ = phat_;
  }

  // resets position and rotation to last saved state
  void reset()  {
    q_ = initialq_;
    phat_ = initialphat_;
  }

  Vec<3> absolutePosOf(Vec<3> relative_pos) {
    return getQ().apply(relative_pos);
  }

  Transformation<> getQ() {
    return q_;
  }

  Vector<double> getPhat()  {
    return phat_;
  }
}; 

class RBS_FEM{

  // references to all bodies added to the system
  std::vector<RigidBody_FEM> bodies_;
  // references to all springs added to the system
  std::vector<Spring> springs_;
   // references to all beams added to the system
  std::vector<Beam> beams_;
   // initialize gravity x, y, z
  Vector<double> gravity_ = {0, 0, 9.81};
  // orderig: first order, second order, first order, second order, ...
  std::vector<double> lag_params;

  Matrix<double> global_mass_inv_{1, 1};


  public:


  //  start
  //  accessor and helper functions for data handling
  //
  std::vector<RigidBody_FEM>& Bodies()  {
    return bodies_;
  }
  RigidBody_FEM& Bodies(size_t i)  {
    return bodies_[i];
  }
  std::vector<Spring>& Springs()  {
    return springs_;
  }
  Spring& Springs(size_t i) {
    return springs_[i];
  }
  std::vector<Beam>& Beams() {
    return beams_;
  }
  Beam& Beams(size_t i)  {
    return beams_[i];
  }
  Vector<double> & Gravity()  {
    return gravity_;
  }
  size_t NumBodies() const  {
    return bodies_.size();
  }
  size_t NumSprings() const {
    return springs_.size();
  }
  size_t NumBeams() const {
    return beams_.size();
  }
  MatrixView<double> getGlobalMassInverse() const {
    return global_mass_inv_; 
  }

  void getState(VectorView<double> out){
    for(int i=0; i<NumBodies(); i++){
      out.Range(i*dim_per_state,i*dim_per_state+12)=bodies_[i].q();
      out.Range(i*dim_per_state + 12, (i+1)*dim_per_state)=bodies_[i].phat();
      std::cout << "getState: " << bodies_[i].phat() << std::endl;
    }
    for (int i=0; i < 2*NumBeams(); i++){
      out(NumBodies()*dim_per_state + i) = lag_params[i];
    }
  }

  void setState(VectorView<double> x)  {
    for(int i=0; i<NumBodies(); i++)  {

      bodies_[i].q() = x.Range(i * dim_per_body, i * dim_per_body + 12);
      bodies_[i].phat() = x.Range(i * dim_per_body + 12, (i + 1) * dim_per_body + 18);
    }

    for (int i=0; i < 2 * NumBeams(); i++)  {

      lag_params[i] = x(NumBodies() * dim_per_body + i);

    }
  }

  void recalcMassMatrixInverse()  {
    Matrix<double> inv(NumBodies()*6, NumBodies()*6);

    for (size_t b=0; b < NumBodies(); b++)  {

      inv.Cols(b*6, 6).Rows(b*6, 6) = inverse(Bodies()[b].Mass_matrix());

    }

    global_mass_inv_ = inv;

  }


  Vector<double> stateToX (VectorView<double> state){
    if((state.Size() - 2*NumBeams())%dim_per_state)
      throw std::invalid_argument("Vector must be in Equation format");

    Vector<double> res(NumBodies() * dim_per_body + 2*NumBeams());
    for(size_t i = 0;i<NumBodies(); i++){
      
      //  copy coordinates
      res.Range(i*dim_per_body,i*dim_per_body+12) = state.Range(i*dim_per_state,i*dim_per_state+12);


      //  coppy impulse
      res.Range(i*dim_per_body+18,i*dim_per_body+24)=state.Range(i*dim_per_state+12,i*dim_per_state+18);
    }

    //  copy Lagrange Parameters
    res.Range(NumBodies()*eq_per_body, res.Size()) = state.Range(NumBodies()*dim_per_state, state.Size());

    return res;
  }


  Vector<double> xToState(VectorView<double> x){
    if((x.Size() - 2*NumBeams()) % eq_per_body)
      throw std::invalid_argument("Vector must be in equation format");

    // the Vector in output format:
    Vector<double> res(NumBodies() * dim_per_state + 2*NumBeams());

    for(size_t i=0; i < NumBodies(); i++){

      //  copy coordinates
      res.Range(i * dim_per_state, i * dim_per_state + 12) = x.Range(i * eq_per_body, i * eq_per_body + 12);
 
      //  copy momentum
      res.Range(i * dim_per_state + 12, i * dim_per_state + 18) = x.Range(i * eq_per_body + 18,i * eq_per_body + 24);
    }

    //  copy Lagrange Parameters
    res.Range(NumBodies() * dim_per_state, res.Size()) = x.Range(NumBodies() * eq_per_body, x.Size());

    return res;
  }


  void saveState()  {
    for (auto& rb: bodies_){
      rb.saveState();
    }
  }

  void reset(){
    for(auto& body: Bodies())
      body.reset();
  }
  //
  //  end
  //  accessor and helper functions for data handling
  

  // start
  // functions for defining and builing the user defined system
  //
  Connector add(RigidBody_FEM& b){

    b.Index() = NumBodies();
    bodies_.push_back(b);

    recalcMassMatrixInverse();
    // b.recalcMassMatrix_large_inverse();
    

    // std::cout << bodies_.size()-1 << std::endl;

    return Connector{ConnectorType::mass, Vector<double>(3), NumBodies()-1};
  }

  // Sets a connector given (x, y, z) relative to body (not absolute), given a body_index 
  Connector SetConnector(double x, double y, double z, size_t body_index) {
    if (body_index >= NumBodies())  {

      throw invalid_argument("Body index out of range");

    } else {

      return Connector{ConnectorType::mass, Vector<double>{x, y, z}, body_index};

    }
  }

  void add(Spring &s){
    //  std::cout << "add spring" << std::endl;
    s.Index() = NumSprings();
    if (s.Connector_a().Type() != ConnectorType::fix)  {
      bodies_[s.Body_index_a()].Springs().push_back(s.Index());
    }
    if (s.Connector_b().Type() != ConnectorType::fix)  {
      bodies_[s.Body_index_b()].Springs().push_back(s.Index());
    }

    //bodies_[s.Body_index_a()].Springs().push_back(s.Index());

    /* if (bodies_[s.Body_index_b()].NumSprings() == 1)  {
      bodies_[s.Body_index_b()].Springs()[0] = s.Index();
    } else  {
      bodies_[s.Body_index_b()].Springs().push_back(s.Index());
    } */

    //bodies_[s.Body_index_b()].Springs().push_back(s.Index());
    
    springs_.push_back(s);
  }

  void add(Beam &b){
    Transformation<> trafo_a = Bodies()[b.Body_index_a()].getQ();
    Transformation<> trafo_b = Bodies()[b.Body_index_b()].getQ();
    b.Length() = Norm(b.Connector_a().absPos(trafo_a)
                    - b.Connector_b().absPos(trafo_b));
    b.Index() = NumBeams();
    //  Notifiy body that they are connected to beam
    if (b.Connector_a().Type() != ConnectorType::fix)  {
      bodies_[b.Body_index_a()].Beams().push_back(b.Index());
    }
    if (b.Connector_b().Type() != ConnectorType::fix)  {
      bodies_[b.Body_index_b()].Beams().push_back(b.Index());
    }
    

    beams_.push_back(b);
    lag_params.push_back(1);
    lag_params.push_back(1);
  }

  Connector addFix(){
    return Connector{ConnectorType::fix, {0,0,0}, 0};
  }

  //  force calculation
  template <typename T>
  Vector<T> gravitation_force(VectorView<T> q, size_t body_index) const {
    Vector<AutoDiffDiff<dim_per_transform, T>> q_diff(dim_per_transform);
    

    for(size_t i = 0; i < dim_per_transform; i++) {
      q_diff(i) = q_diff(i);
      q_diff(i).DValue(i) = 1;
    }

    //  std::cout << q << std::endl; 

    auto t = Transformation<AutoDiffDiff<dim_per_transform, T>>(q_diff);
    VectorView<T> res(dim_per_transform, ((-1)*bodies_[body_index].mass() * t.apply(bodies_[body_index].center()) * gravity_).DValue());
    //std::cout << ((-1)*bodies_[body_index].mass() * t.apply(bodies_[body_index].center()) * gravity_).Value() << std::endl;
    return res;
  }

  template<typename T, typename S>
  void force(VectorView<T>& x, VectorView<S>& f, size_t body_index) const {
    f = 0;
    f.Range(0, dim_per_transform) -= gravitation_force(x.Range(dim_per_body * body_index, 
                                        dim_per_body * body_index + dim_per_transform), body_index);

    for (size_t i: bodies_[body_index].Springs()) {
      Spring spr = springs_[i];

      size_t body_index_b = (body_index != spr.Body_index_a())? spr.Body_index_a() : spr.Body_index_b();
      body_index_b = (spr.Connector_a().Type() == ConnectorType::fix)? spr.Body_index_b() : body_index_b;
      
      Vector<T> s_f = spr.force(x.Range(dim_per_body * body_index, 
                                        dim_per_body * body_index + dim_per_transform),
                                x.Range(dim_per_body * body_index_b,
                                        dim_per_body * body_index_b + dim_per_transform));

      //  Adding force for body a
      f.Range(0, dim_per_transform) -= s_f.Range(0, dim_per_transform);
    }

    for (size_t i: bodies_[body_index].Beams()) {
      Beam bm = beams_[i];

      size_t body_index_b = (body_index != bm.Body_index_a())? bm.Body_index_a() : bm.Body_index_b();
      body_index_b = (bm.Connector_a().Type() == ConnectorType::fix)? bm.Body_index_b() : body_index_b;
      
      Vector<T> b_f = bm.force(x.Range(dim_per_body * body_index, 
                                        dim_per_body * body_index + dim_per_transform),
                                x.Range(dim_per_body * body_index_b,
                                        dim_per_body * body_index_b + dim_per_transform),
                                x(dim_per_body*NumBodies() + 2 * i));

      //  Adding force for body a
      f.Range(0, dim_per_transform) -= b_f.Range(0, dim_per_transform);
    }

    //  std::cout << "Body Index: " << body_index << "Force: " << f << std::endl;
  }

  // first beam constraint
  template<typename T> 
  T g(VectorView<T>& x, size_t beam_index)  {
    Beam bm = beams_[beam_index];
    // calclation of potential and force
    size_t k = bm.Connector_a().Body_index();
    size_t l = bm.Connector_b().Body_index();

    Vec<3, T> pos1 = bm.Connector_a().absPos(x.Range(dim_per_body * k, dim_per_body * k + dim_per_transform));
    Vec<3, T> pos2 = bm.Connector_b().absPos(x.Range(dim_per_body * l, dim_per_body * l + dim_per_transform));

    T norm = Norm(pos1-pos2) - bm.Length();
    return (norm * norm);
  }

  //  derivative of first constraint
  template<typename T>
  Vector<T> G(VectorView<T> q_a, VectorView<T> q_b, Beam& bm) {
    
    Vector<AutoDiffDiff<2*dim_per_transform, T>> res(1);

    Vector<AutoDiffDiff<2*dim_per_transform, T>> x_diff(2*dim_per_transform);

    for (size_t i = 0; i < dim_per_transform; i++)  {
      // first dim_per_transform are for body a
      // setting Values and Differential Index
      x_diff(i) = q_a(i);
      x_diff(i).DValue(i) = 1;
      x_diff(dim_per_transform + i) = q_b(i);
      x_diff(dim_per_transform + i).DValue(dim_per_transform + i) = 1;
    }

    // calclation of potential and force
    Vec<3, AutoDiffDiff<2*dim_per_transform, T>> pos1 = bm.Connector_a().absPos(x_diff.Range(0, dim_per_transform));
    Vec<3, AutoDiffDiff<2*dim_per_transform, T>> pos2 = bm.Connector_b().absPos(x_diff.Range(dim_per_transform, 2*dim_per_transform));

    AutoDiffDiff<2*dim_per_transform, T> norm = Norm(pos1-pos2) - bm.Length();
    res(0) = (norm * norm);
    
    return VectorView(2*dim_per_transform, res(0).DValue());
  }

  //  velocityconstraint for one beam
  template<typename T>
  T velocity_constraint(VectorView<T>& x, size_t beam_index) {

    Beam bm = beams_[beam_index];

    size_t body_index_b = bm.Body_index_b();
    size_t body_index_a = bm.Body_index_a();

    Vector<T> G_i = G(x.Range(body_index_a * dim_per_body, body_index_a * dim_per_body + dim_per_transform), 
                      x.Range(dim_per_body * body_index_b, dim_per_body * body_index_b + dim_per_transform), bm);

    Vector<T> temp(2*dim_per_transform);

    temp.Range(0, 3) = 1./bodies_[body_index_a].mass() * x.Range(body_index_a * dim_per_body + 18, body_index_a * dim_per_body + 21);

    Matrix<T> B_a = bodies_[body_index_a].Mass_matrix_large_invers()*AsMatrix(x.Range(body_index_a * dim_per_body + 3, 
                                                                                                body_index_a * dim_per_body + dim_per_transform), 3, 3)
                                 *hat_map(x.Range(body_index_a * dim_per_body + 21, body_index_a * dim_per_body + 24));
    
    Matrix<T> B_b = bodies_[body_index_b].Mass_matrix_large_invers()*AsMatrix(x.Range(body_index_b * dim_per_body + 3, 
                                                                                                body_index_b * dim_per_body + dim_per_transform), 3, 3)
                                 *hat_map(x.Range(body_index_b * dim_per_body + 21, body_index_b * dim_per_body + 24));

    temp.Range(3, 12) = AsVector(B_a);
      
    temp.Range(12, 15) = 1./bodies_[body_index_b].mass() * x.Range(body_index_b*dim_per_body + 18, body_index_b*dim_per_body + 21);

    temp.Range(15, 24) = AsVector(B_b);

    return G_i*temp;
  }

  //  derivative of all velocity constraints with respect to one body
  template<typename T, typename S>
  void dvelocity_constraint(VectorView<T>& x, VectorView<S>& f, size_t body_index)  {

    AutoDiffDiff<dim_per_state, T> res;
    Vector<AutoDiffDiff<dim_per_state, T>> x_diff(dim_per_state);

    for (size_t i = 0; i <dim_per_transform; i++) {
      x_diff(i) = x(body_index*dim_per_body + i);
      x_diff(i).DValue(i) = 1;
    }

    for (size_t i = 0; i < 6; i++)  {
      x_diff(dim_per_transform + i) = x(body_index*dim_per_body + dim_per_transform + i);
      x_diff(dim_per_transform + i) = 1;
    }


    for (size_t i: bodies_[body_index].Beams()) {
      Beam bm = beams_[i];

      size_t body_index_b = (body_index != bm.Body_index_a())? bm.Body_index_a() : bm.Body_index_b();
      body_index_b = (bm.Connector_a().Type() == ConnectorType::fix)? bm.Body_index_b() : body_index_b;

      Vector<AutoDiffDiff<dim_per_state, T>> q_b_diff = x.Range(dim_per_body * body_index_b,
                                                                dim_per_body * body_index_b + dim_per_transform);

      Vector<AutoDiffDiff<dim_per_state, T>> G_i = G(x_diff.Range(0, dim_per_transform), q_b_diff, bm);

      Vector<AutoDiffDiff<dim_per_state, T>> temp(2*dim_per_transform);

      temp.Range(0, 3) = 1./bodies_[body_index].mass() * x_diff.Range(dim_per_transform, dim_per_transform + 3);

      Matrix<AutoDiffDiff<dim_per_state, T>> B_a = bodies_[body_index].Mass_matrix_large_invers()*AsMatrix(x_diff.Range(3, dim_per_transform), 3, 3)
                          *hat_map(x_diff.Range(dim_per_transform + 3, 2 * dim_per_transform));

      Matrix<AutoDiffDiff<dim_per_state, T>> B_b = bodies_[body_index_b].Mass_matrix_large_invers()*AsMatrix(x.Range(body_index_b*dim_per_body + 3, 
                                                                                                body_index_b*dim_per_body + dim_per_transform), 3, 3)
                          *hat_map(x.Range(body_index_b*dim_per_body + 21, body_index_b*dim_per_body + 24));

      temp.Range(3, 12) = AsVector(B_a);

      
      
      temp.Range(12, 15) = 1./bodies_[body_index_b].mass() * x.Range(body_index_b*dim_per_body + 18, body_index_b*dim_per_body + 21);

      temp.Range(15, 24) = AsVector(B_b);


      //std::cout << temp << std::endl;
      
      AutoDiffDiff<dim_per_state, T> rs = G_i*temp;

      //std::cout << "res: " << rs << std::endl;

      res = res + x(NumBodies()*dim_per_body + 2*bm.Index() + 1)*rs;
      //std::cout << res << std::endl;

    }
  

    for (size_t i = 0; i < dim_per_state; i++) {
      f(i) = res.DValue(i);
    }
  }

};

#endif