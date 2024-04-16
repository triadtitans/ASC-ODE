#include <nonlinfunc.h>
#include <math.h>
#include <matrix.h>
#include <vector.h>
#include <ode.h>
#include "../src/autodiffdiff.h"

constexpr size_t dim_per_body = 18;

class RhsRigidBody;

using namespace ASC_bla;
using namespace ASC_ode;

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

  Vector<T> getTranslation(){
    return Vector<T>{q_(0),q_(4),q_(8)};
  }

  Matrix<T> getRotation() {
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


class RigidBody {
  Vector<double> q_;
  Vector<double> dq_;
  Vector<double> phat_;
  Vector<double> initialq_;
  Vector<double> initialdq_;
  Vector<double> initialphat_;

  Matrix<double> T;

  Vec<3> center_of_mass_={0,0,0};
  Matrix<double> inertia_;
  double mass_=1;
  std::shared_ptr<NonlinearFunction> mass_function;
public:

  RigidBody(Vector<double> q,Vector<double> dq,Vector<double> phat,double mass, Vec<3> center_of_mass,MatrixExpr<T>& inertia)
        : q_(q),dq_(dq), phat_(phat), initialq_(q), initialdq_(dq), initialphat_(phat), center_of_mass_(center_of_mass),
          mass_(mass), inertia_(inertia), T(6, 12) {
    if(inertia.Width() != 3) throw std::invalid_argument("Inertia matrix must be 3x3");
    if(inertia.Height() != 3) throw std::invalid_argument("Inertia matrix must be 3x3");
    if(q.Size() != dim_per_body) throw std::invalid_argument("q Vector must match mass matrix");
    if(dq.Size() != dim_per_body) throw std::invalid_argument("q Vector must match mass matrix");
    if(phat.Size() !=6) throw std::invalid_argument("q Vector must match mass matrix");

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
  void setDq(Transformation<> t){dq_=t.q_;}
  void setPhat(Vector<double> v){phat_=T*v;}

  void setMass(Matrix<double> m){mass_function=std::make_shared<LinearFunction>(m);}
  std::shared_ptr<NonlinearFunction> getMassFunc(){return mass_function;}
  // saves a state for the reset button
  void saveState(){
    initialq_ = q_;
    initialdq_ = dq_;
    initialddq_ = ddq_;
  }
  // resets position and rotation to last saved state
  void reset(){
    q_ = initialq_;
    dq_ = initialdq_;
    ddq_ = initialddq_;
  }

  Vec<3> absolutePosOf(Vec<3> relative_pos){
    return getQ().apply(relative_pos);
  }

  Transformation<> getQ(){return q_;}
  Transformation<> getDq(){return dq_;}
  Transformation<> getDdq(){return ddq_;}

  void simulate(double tend, double steps, std::function<void(double,VectorView<double>)> callback = nullptr ){
    /* std::shared_ptr<RhsRigidBody> rhs = std::make_shared<RhsRigidBody>(*this);
    //std::shared_ptr<NumericDerivative> dlagrange = std::make_shared<NumericDerivative>(rhs);
    std::shared_ptr<Derivative> dlagrange = std::make_shared<Derivative>(rhs);
   
    SolveODE_Alpha (tend, steps, 0.8, q_, dq_, ddq_, dlagrange, mass_function, callback); */


  } 
};

class EQRigidBody : public NonlinearFunction
{
  double h_;
  RigidBody b_;
  public:
  //Just a reminder
  EQRigidBody(RigidBody b, double h) h_(h), b_(b) {};
  size_t DimX() const  { return 39; } //???
  size_t DimF() const  { return 48; } //???

  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    // x (a(3), B(row_maj)(9), v_trans(3), v_skew(3), phat(6), p(6), Bhalf(9))
    Vector<double> anew = x.Range(0, 3);
    MatrixView<double> Bnew = asMatrix(x.Range(3, 12), 3, 3);
    Vector<double> vtrans = x.Range(12, 15);
    Vector<double> vskew = x.Range(15, 18);
    Vector<double> phat = x.Range(18, 24);
    Vector<double> p = x.Range(24, 30);
    MatrixView<double> Bhalf = asMatrix(x.Range(30, 39), 3, 3);

    Vector<double> aold = b_.getQ().getTranslation();
    Matrix<double> Bold = b_.getQ().getRotation();

    f.Range(0, 3) = (1/h_)*(anew - aold) - vtrans;
    MatrixView<double> f_temp(asMatrix(f.Range(3, 12)));
    ftemp = (1/h_)*(Bnew - Bold) - Bhalf*vskew; // Bold oder B_{1/2}???
    // M*(vtrans, vskew) = phatold
    MatrixView<double> ftemp(asMatrix(f.Range(12, 21)));
    ftemp = Bold*phatold - Bhalf*pold; // pold oder phatold???
    MatrixView<double> ftemp(asMatrix(f.Range(21, 30)));
    ftemp = Bnew*phat - Bhalf*pold;
    MatrixView<double> ftemp(asMatrix(f.Range(30, 39)));
    ftemp = Transpose(Bnew)*Bnew - Diagonal(3, 1);
    MatrixView<double> ftemp(asMatrix(f.Range(39, 48)));
    ftemp = Transpose(Bhalf)*Bhalf - Diagonal(3, 1);
  }

  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    
  }
};

  



/* 
inline double mass_matrix_data[18*18] ={
1.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
-0.0, 0.08333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
-0.0, 0.0, 0.08333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.08333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, 1.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, -0.0, 0.08333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, -0.0, 0.0, 0.08333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.08333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, 0.08333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, 0.0, 0.08333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.08333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

inline double inertia_matrix_data[3*3] ={
0.08333, 0.0, 0.0, \
 0.0, 0.08333, 0.0,\
 0.0, 0.0, 0.08333, \
}; */

