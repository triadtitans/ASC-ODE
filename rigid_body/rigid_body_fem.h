#include <nonlinfunc.h>
#include <math.h>
#include <matrix.h>
#include <vector.h>
#include <ode.h>
#include "../src/autodiffdiff.h"

constexpr size_t dim_per_body = 18;

class RigidBody;

using namespace ASC_bla;
using namespace ASC_ode;

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
  Vector<T> q_;

  Transformation(Vector<T> q):q_(q){}
  Transformation():q_(18){}

  Vec<3, T> apply(Vec<3, T> pos){
    Vec<3, T> trans{q_(0),q_(4),q_(8)};
    Matrix<T> rot = AsMatrix(q_,3,4).Cols(1,3);
    return trans + rot*pos;
  }

  void setTranslation(T a, T b, T c){q_(0)=a;q_(4)=b;q_(8)=c;}
  void setRotation(int i, int j, double r){
    if(i>2||j>2 || j<0 || i <0) throw std::invalid_argument("Rotation Matrix is 3x3");
    q_((4*i + 1) + j)=r; // note the ordering of q
  }
  template<typename S>
  void setRotation(const MatrixExpr<S>& B){
    for (size_t i=0; i < 3; i++){
      for (size_t j=0; j < 3; j++){
        setRotation(i, j, B(i, j));
      }
    }
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


class EQRigidBody : public NonlinearFunction
{
  double h_;
  Transformation<double> Q_;
  Vector<double> phatold;
  public:
  EQRigidBody(Transformation<double> Q, Vector<double> Phat, double h): h_(h), Q_(Q), phatold(Phat) {};
  size_t DimX() const  { return 30; }
  size_t DimF() const  { return 30; }

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

    // I
    f.Range(0, 3) = (1/h_)*(anew - aold) - vtrans;
    f.Range(3, 6) = (1/h_)*inv_hat_map(Transpose(Bhalf)*(Bnew - Bold)) - vtrans;

    // II
    Vector<double> vhat(6);
    vhat.Range(0, 3) = vtrans;
    vhat.Range(3, 6) = vskew;

    // for the moment, M shall be id
    f.Range(6, 12) = vhat - p;

    // III
    
    f.Range(12, 15) = (2/h_)*(phat.Range(0, 3) - p.Range(0, 3));
    f.Range(15, 18) = inv_hat_map((2/h_)*(Transpose(Bnew) * (Bnew*hat_map(phat.Range(3, 6)) - Bhalf*hat_map(p.Range(3, 6)))));
    f.Range(18, 21) = (2/h_)*(p.Range(0, 3) - phatold.Range(0, 3));
    f.Range(21, 24) = inv_hat_map((2/h_)*(Transpose(Bold) * (Bhalf*hat_map(p.Range(3, 6)) - Bold*hat_map(phatold.Range(3, 6)))));

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


class RigidBody {
  Vector<double> q_;
  Vector<double> phat_;
  Vector<double> initialq_;
  Vector<double> initialphat_;

  Matrix<double> P_;

  Vec<3> center_of_mass_={0,0,0};
  Matrix<double> inertia_;
  double mass_=1;
  std::shared_ptr<NonlinearFunction> mass_function;
public:

  template<typename T>
  RigidBody(Vector<double> q,Vector<double> phat,double mass, Vec<3> center_of_mass,MatrixExpr<T>& inertia)
        : q_(q), phat_(phat), initialq_(q), initialphat_(phat), center_of_mass_(center_of_mass),
          mass_(mass), inertia_(inertia), P_(6, 12) {
    if(inertia.Width() != 3) throw std::invalid_argument("Inertia matrix must be 3x3");
    if(inertia.Height() != 3) throw std::invalid_argument("Inertia matrix must be 3x3");
    if(q.Size() != dim_per_body) throw std::invalid_argument("q Vector must match mass matrix");
    if(phat.Size() != 6) throw std::invalid_argument("q Vector must match mass matrix");

  }

  /* 
  RigidBody()
        :   mass_function(std::make_shared<LinearFunction>(Matrix(18,18))),
          q_(18),dq_(18),ddq_(18), initialq_(18),initialdq_(18),initialddq_(18),inertia_(3,3),center_of_mass_{0,0,0}{
    q_(1)=1;q_(6)=1;q_(11)=1;
    recalcMassMatrix();
  }
 */
  double& mass(){return mass_;}
  Vec<3>& center(){return center_of_mass_;}
  Matrix<double>& inertia(){return inertia_;}

  /* 
  void recalcMassMatrix(){
    mass_function = nullptr;
    //Times 2 because of derivative of x * Ax in x is Ax + A(T)x
    auto diag_function=std::make_shared<LinearFunction>(2*diagonal_block_from_inertia(inertia_,center_of_mass_,mass_));
    auto block_func =std::make_shared<BlockFunction>(diag_function,3);
    auto mass = std::make_shared<StackedFunction>();
    mass->addFunction(block_func);
    Vector<double> zero(6);
    std::shared_ptr<NonlinearFunction> const_zero = std::make_shared<ConstantFunction>(zero);
    mass->addFunction(const_zero);

    mass_function = mass;
  }
  */

  void setQ(Transformation<> t){q_=t.q_;}
  void setPhat(Vector<double> v){phat_=v;}

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
    return getQ().apply(relative_pos);
  }

  Transformation<> getQ(){return q_;}
  Vector<double> getPhat(){return phat_;}

  void simulate(double tend, double steps, std::function<void(int,double,VectorView<double>)> callback = nullptr ){
    Vector<double> state(30);
    state.Range(0,3)=getQ().getTranslation();
    AsMatrix(state.Range(3, 12), 3, 3) = getQ().getRotation();
    state.Range(18, 24) = phat_;
    for (size_t step=0; step < steps; step++){
      auto q = getQ();
      auto p = getPhat();
      std::shared_ptr<EQRigidBody> eq = std::make_shared<EQRigidBody>(q, p, tend/steps);
      NewtonSolver(eq, state, 1e-10, 10, callback);
      Transformation<double> Q = getQ();
      Q.setTranslation(state(0),state(1),state(2));
      Q.setRotation(AsMatrix(state.Range(3, 12), 3, 3));
      setQ(Q);
      phat_ = state.Range(18, 24);
    }
  } 
};
