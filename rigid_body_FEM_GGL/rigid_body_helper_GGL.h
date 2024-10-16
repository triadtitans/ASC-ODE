#ifndef rb_helper_GGL
#define rb_helper_GGL

#include <nonlinfunc.h>
#include <math.h>
#include <matrix.h>
#include <vector.h>
#include <ode.h>
#include "../src/autodiffdiff.h"
#include <cmath>
#include <chrono>


#ifdef PYBIND11_MODULE
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif


constexpr size_t dim_per_transform = 12;
constexpr size_t dim_per_body = 30;
constexpr size_t dim_per_state = 18;
constexpr size_t eq_per_body = 30;

using namespace ASC_bla;
using namespace ASC_ode;


// convert vector v in R³ to skew-symmetric matrix
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

// convert vector v in R³ to vector h = (v[0,1,2], hat_map(v[3,4,5]))
template<typename T>
auto hat_map_vector(const VectorView<T>& v){
  Vector<T> h(dim_per_transform);
  h.Range(0, 3) = v.Range(0, 3);
  h.Range(3, dim_per_transform) = AsVector(hat_map(v.Range(3, dim_per_transform)));
  return h;
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
  // ordering of q_: (trans_x, trans_y, trans_z, rot_xx, rot_xy, rot_xz, rot_yx, rot_yy, rot_yz, rot_zx, rot_zy, rot_zz)
  Vector<T> q_;

  Transformation(Vector<T> q):q_(q){}
  Transformation():q_(12){
    q_(3)=1;
    q_(7)=1;
    q_(11)=1;
  }

  // applys transformation to given position
  Vec<3, T> apply(Vec<3, T> pos){
    Vec<3, T> trans{q_(0),q_(1),q_(2)};
    Matrix<T> rot = AsMatrix(q_.Range(3, 12),3,3);
    return trans + rot*pos;
  }

  // set Translation with single entries
  void setTranslation(T a, T b, T c){q_(0)=a;q_(1)=b;q_(2)=c;}

  // set Translation with Vector
  void setTranslation(VectorView<T> v) {q_(0)=v(0), q_(1)=v(1); q_(2)=v(2);}

  // set Rotation entry-wise
  void setRotation(size_t i, size_t j, T r){
    if(i>2||j>2 || j<0 || i <0) throw std::invalid_argument("Rotation Matrix is 3x3");
    q_((3+ 3*i) + j)=r; 
  }


  //set Rotation from Matrix
  void setRotation(MatrixView<T> B){
    for (size_t i=0; i < 3; i++){
      for (size_t j=0; j < 3; j++){
        q_((3+ 3*i) + j) = B(i, j);
      }
    }
  }

  // set rotation matrix from x, y or z degrees
  void setRotationDeg(int axis, T deg){
    Matrix<T> R = makeRotationMatrix3<T>(axis, deg);
    setRotation(R);
  }

  Vector<T> getTranslation() const{
    Vector<T> v(3);
    v(0)=q_(0);
    v(1)=q_(1);
    v(2)=q_(2);
    return v;
  }
  Matrix<T> getRotation() const {
    Matrix<T> B(3, 3);
    for (size_t i=0; i < 3; i++){
      for (size_t j=0; j < 3; j++){
        B(i, j) = q_(3 + 3*i + j);
      }
    }
    return B;
  }
  double getAngle(){
    Matrix<double> rot = getRotation();
    double tr = rot(0, 0) + rot(1, 1) + rot(2, 2);
    return radToDeg(std::acos((tr - 1)/2));
  }
};

std::ostream& operator<<(std::ostream& oss, const Transformation<double>& t){
    oss<<std::fixed << " Translation: \t" << t.q_(0) << " ," << t.q_(1) << ", "<<", " << t.q_(2) << std::endl
                    << " Rotation: \t" << t.q_(3) << " ," << t.q_(4) << ", "<<", " << t.q_(5) <<  std::endl
                    << " \t\t" << t.q_(6) << " ," << t.q_(7) << ", "<<", " << t.q_(8) <<  std::endl
                    << "\t\t" << t.q_(9) << " ," << t.q_(10) << ", "<<", " << t.q_(11) << std::endl; 
  return oss;
}

// SE(3)
template<typename T>
Transformation<T> operator*(Transformation<T> t1, Transformation<T> t2){
  Transformation<T> res;
  Vector<T> trans = t1.getTranslation() + t1.getRotation()*t2.getTranslation();
  res.setTranslation(trans(0), trans(1), trans(2));
  res.setRotation(t1.getRotation() * t2.getRotation());

  return res;
}

// ℝ x SO(3)
template<typename T>
Transformation<T> operator+(Transformation<T> t1, Transformation<T> t2){
  Transformation<T> res;
  Vector<T> trans = t1.getTranslation() + t1.getRotation()*t2.getTranslation();
  res.setTranslation(trans(0), trans(1), trans(2));
  res.setRotation(t1.getRotation() * t2.getRotation());

  return res;
}



enum class ConnectorType {mass, fix};

// manages connections to bodies with Springs etc.
class Connector{
  ConnectorType t_;
  size_t body_index_; // index of body in system's _bodies that the Connector is relative to; irrelevant for fixes

  // size_t num;
  Vec<3> pos_={0,0,0};

  public:

  Connector (ConnectorType t, Vec<3> pos, size_t body_index): pos_(pos), body_index_(body_index), t_(t) {;}

  size_t& Body_index()  {
    return body_index_;
  }

  ConnectorType& Type()  {
    return t_;
  }

  Vec<3>& Pos()  {
    return pos_;
  }

  double& Pos(size_t i)  {
    return pos_(i);
  }

  template<typename T>
  Vec<3, T> absPos(Transformation<T> trafo) {
    if(t_ == ConnectorType::fix){

      Vec<3, T> pos_t(pos_);
      return pos_t;

    }  else  {

      Vec<3, T> res = trafo.apply(pos_);
      return res;

    }
  }

  template<typename T>
  Vec<3, T> absPos(Vector<T> a, Matrix<T> B) const{
    
    if(t_ == ConnectorType::fix){

      Vec<3, T> pos_t(pos_);
      return pos_t;

    } else  {

      Transformation<T> tr;
      tr.setTranslation(a(0), a(1), a(2));
      tr.setRotation(B);
      Vec<3, T> res = tr.apply(pos_);
      return res;

    }
  }

  template<typename T>
  Vec<3, T> absPos(VectorView<T> q) const
  {
    
    if(t_ == ConnectorType::fix) {

      Vec<3, T> pos_t(pos_);
      return pos_t;

    }
    else  {

      Transformation<T> tr(q);
      Vec<3, T> res = tr.apply(pos_);
      return res;

    }
  }
};

std::ostream& operator<<(std::ostream& oss, Connector& c){
  if (c.Type() == ConnectorType::fix) {
    oss<< std::fixed << "Type: fix" << std::endl;
  } else  {
    oss<< std::fixed << "Type: mass" << std::endl;
  }
  oss<< std::fixed << "Body_index: " << c.Body_index() << std::endl
                   << "Position: " << c.Pos() << std::endl;
  return oss;
}


// Handling position and force for Springs added to the system
class Spring{
  double length_;
  double stiffness_;
  size_t index_;
  Connector a_;
  Connector b_;

  public:

  Spring(double length, double stiffness, Connector a, Connector b): length_(length), stiffness_(stiffness), a_(a), b_(b) {};

  size_t& Index()  {
    return index_;
  }

  size_t Body_index_a() {
    return a_.Body_index();
  }

  size_t Body_index_b() {
    return b_.Body_index();
  }

  double& Length()  {
    return length_;
  }

  double& Stiffness()  {
    return stiffness_;
  }

  Connector& Connector_a() {
    return a_;
  }

  Connector& Connector_b() {
    return b_;
  }

  template<typename T>
  Vector<T> force(VectorView<T> q_a, VectorView<T> q_b, bool diff_index) {
    // bool diff_index = True -> differentiate after q_a
    // else differentiate after q_b
    
    Vector<AutoDiffDiff<dim_per_transform, T>> res(1);

    Vector<AutoDiffDiff<dim_per_transform, T>> x_diff(dim_per_transform);

    if (diff_index) {
      for (size_t i = 0; i < dim_per_transform; i++)  {
        // first dim_per_transform are for body a
        // setting Values and Differential Index
        x_diff(i) = q_a(i);
        x_diff(i).DValue(i) = 1;
      }

      Vec<3, AutoDiffDiff<dim_per_transform, T>> pos1 = Connector_a().absPos(x_diff.Range(0, dim_per_transform));
      Vec<3, T> pos2 = Connector_b().absPos(q_b);

      AutoDiffDiff<dim_per_transform, T> norm = Norm(pos1-pos2) - length_;
      res(0) = (1/2.0)*Stiffness()*(norm * norm);
      Vector<T> res_f = res(0).DValue_vec();
      //  std::cout << "Sprint_force: " << res_f << std::endl;
      return res_f;
      
    } else  {
      for (size_t i = 0; i < dim_per_transform; i++)  {
        // first dim_per_transform are for body a
        // setting Values and Differential Index
        x_diff(i) = q_b(i);
        x_diff(i).DValue(i) = 1;
      }

      Vec<3, T> pos1 = Connector_a().absPos(q_a);
      Vec<3, AutoDiffDiff<dim_per_transform, T>> pos2 = Connector_b().absPos(x_diff.Range(0, dim_per_transform));

      AutoDiffDiff<dim_per_transform, T> norm = Norm(pos1-pos2) - length_;
      res(0) = (1/2.0)*Stiffness()*(norm * norm);
      Vector<T> res_f = res(0).DValue_vec();
      //  std::cout << "Sprint_force: " << res_f << std::endl;
      return res_f;
      
    }
  }
};

std::ostream& operator<<(std::ostream& oss, Spring& sp){
  oss<<std::fixed << "Index: " << sp.Index() << std::endl
                  << "Length: " << sp.Length() << std::endl
                  << "Stiffness: " << sp.Stiffness() << std::endl
                  << "Connector A: " << std::endl << sp.Connector_a() <<  std::endl
                  << "Connector B: " << std::endl << sp.Connector_b() <<  std::endl;
  return oss;
}


// Handling position and forces for Beams added to the system
class Beam  {
  double length_;
  size_t index_;
  Connector a_;
  Connector b_;

  public:

  Beam(Connector a, Connector b): index_(0), a_(a), b_(b)  {}

  size_t& Index()  {
    return index_;
  }

  size_t Body_index_a() {
    return a_.Body_index();
  }

  size_t  Body_index_b() {
    return b_.Body_index();
  }

  Connector& Connector_a() {
    return a_;
  }

  Connector& Connector_b() {
    return b_;
  }

  double& Length() {
    return length_;
  }


  //  takes a vector with oredring (values_body_1, values_body_2, ... , values_body_n, 
  //  lambda_1_body_1, lambda_2_bdy_2, ... , lambda_1_body_n, lambda_2_body_n)
  //  ordering of body_values are assumed to be (trans, rotation (rowmajor), ...)
  template<typename T>
  Vector<T> force(VectorView<T> q_a, VectorView<T> q_b, T lambda, bool diff_index) {
    // bool diff_index = True -> differentiate after q_a
    // else differentiate after q_b
    
    Vector<AutoDiffDiff<dim_per_transform, T>> res(1);

    Vector<AutoDiffDiff<dim_per_transform, T>> x_diff(dim_per_transform);

    if (diff_index) {
      for (size_t i = 0; i < dim_per_transform; i++)  {
        // first dim_per_transform are for body a
        // setting Values and Differential Index
        x_diff(i) = q_a(i);
        x_diff(i).DValue(i) = 1;
      }

      Vec<3, AutoDiffDiff<dim_per_transform, T>> pos1 = Connector_a().absPos(x_diff.Range(0, dim_per_transform));
      Vec<3, T> pos2 = Connector_b().absPos(q_b);

      res(0) =lambda*((pos1 - pos2)*(pos1 - pos2) - Length()*Length());
      //std::cout << "pos1: " << pos1 << std::endl;
      //std::cout << "pos2: " << pos2 << std::endl;
      //std::cout << "length: " << length_ << std::endl;
      //std::cout << "lambda: " << lambda << std::endl;
    
      return res(0).DValue_vec();
      
    } else  {
      for (size_t i = 0; i < dim_per_transform; i++)  {
        // first dim_per_transform are for body a
        // setting Values and Differential Index
        x_diff(i) = q_b(i);
        x_diff(i).DValue(i) = 1;
      }

      Vec<3, T> pos1 = Connector_a().absPos(q_a);
      Vec<3, AutoDiffDiff<dim_per_transform, T>> pos2 = Connector_b().absPos(x_diff.Range(0, dim_per_transform));
      //AutoDiffDiff<dim_per_transform, T> res = ((pos1-pos2)*(pos1-pos2) - beam.length*beam.length)*((pos1-pos2)*(pos1-pos2) - beam.length*beam.length);
      //AutoDiffDiff<dim_per_transform, T> norm = Norm(pos1-pos2) - Length();
      res(0) = lambda*((pos1 - pos2)*(pos1 - pos2) - Length()*Length());
      //res(0) = lambda*((pos1-pos2)*(pos1-pos2) - length_*length_)*((pos1-pos2)*(pos1-pos2) - length_*length_);
      //std::cout << "pos1: " << pos1 << std::endl;
      //std::cout << "pos2: " << pos2 << std::endl;
      //std::cout << "length: " << length_ << std::endl;
      //std::cout << "lambda: " << lambda << std::endl;
      return res(0).DValue_vec();
      
    }
    
    // lambda parameters
    // setting Values and Differential Index
    //x_diff(dim_per_transform) = lambda;
    //x_diff(dim_per_transform).DValue(dim_per_transform) = 1;
  }
};

std::ostream& operator<<(std::ostream& oss, Beam& bm){
  oss<<std::fixed << "Index: " << bm.Index() << std::endl
                  << "Length: " << bm.Length() << std::endl
                  << "Connector A: " << std::endl << bm.Connector_a() <<  std::endl
                  << "Connector B: " << std::endl << bm.Connector_b() <<  std::endl;
  return oss;
}

#endif