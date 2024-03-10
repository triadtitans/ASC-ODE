#include <nonlinfunc.h>
#include <math.h>
#include <matrix.h>
#include <vector.h>
#include <ode.h>

class RhsRigidBody;

using namespace ASC_bla;
using namespace ASC_ode;

struct Transformation{
  Vector<double> q_;
  Transformation(Vector<double> q):q_(q){};
  Transformation():q_(18){}
  void setTranslation(double a, double b, double c){q_(0)=a;q_(1)=b;q_(2)=c;}
  void setRotation(int i, int j, double r){
    if(i>2||j>2 || j<0 || i <0) throw std::invalid_argument("Rotation Matrix is 3x3");
    q_(3+i*3+j)=r;
  }                               
};

  std::ostream& operator<<(std::ostream& oss, const Transformation& t){
      oss<<std::fixed << " Translation: \t" << t.q_(0) << " ," << t.q_(1) << ", "<<", " << t.q_(2) << std::endl
                      << " Rotation: \t" << t.q_(3) << " ," << t.q_(4) << ", "<<", " << t.q_(5) <<  std::endl
                     << " \t\t" << t.q_(6) << " ," << t.q_(7) << ", "<<", " << t.q_(8) <<  std::endl
                      << "\t\t" << t.q_(9) << " ," << t.q_(10) << ", "<<", " << t.q_(11) << std::endl; 
    return oss;
  }    

struct MassMatrix{
  Matrix<double> m_;
  MassMatrix(Matrix<double> m): m_(m){}
  MassMatrix(): m_(18,18){}
  void set(int i, int j, double r){
    m_(i,j)=r;
  }
};


class RigidBody {
  Vector<double> q_;
  Vector<double> dq_;
  Vector<double> ddq_;
  int dim_;
  std::shared_ptr<LinearFunction> mass_function;
public:
  template <typename T>
  RigidBody(MatrixExpr<T>& m,Vector<double> q,Vector<double> dq,Vector<double> ddq)
        : dim_(m.Height()), mass_function(std::make_shared<LinearFunction>(m)),
          q_(q),dq_(dq),ddq_(ddq){
    if(m.Width() != dim_) throw std::invalid_argument("Mass matrix must be square");
    if(q.Size() != dim_) throw std::invalid_argument("q Vector must match mass matrix");
    if(dq.Size() != dim_) throw std::invalid_argument("q Vector must match mass matrix");
    if(ddq.Size() != dim_) throw std::invalid_argument("q Vector must match mass matrix");
  }

  RigidBody()
        :  dim_(18), mass_function(std::make_shared<LinearFunction>(Matrix(18,18))),
          q_(18),dq_(18),ddq_(18){
    q_(3)=1;q_(7)=1;q_(11)=1;
  }

  void setQ(Transformation t){q_=t.q_;}
  void setDq(Transformation t){dq_=t.q_;}
  void setDdq(Transformation t){ddq_=t.q_;}

  void setMass(MassMatrix m){mass_function=std::make_shared<LinearFunction>(m.m_);}

  Transformation getQ(){return q_;}
  Transformation getDq(){return dq_;}
  Transformation getDdq(){return ddq_;}

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

inline double mass_matrix_data[18*18] = 
{0.5, 0., 0., 0.25, 0.25, 0.25, 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0.5, 0., 0., 0., 0., 0.25, 0.25, 0.25, 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0., 0., 0., 0., 0., 0., 0.25, \
0.25, 0.25, 0., 0., 0., 0., 0., 0., 0.25, 0., 0., 0.166667, 0.125, \
0.125, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0., 0., \
0.125, 0.166667, 0.125, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0.25, 0., 0., 0.125, 0.125, 0.166667, 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0.25, 0., 0., 0., 0., 0.166667, 0.125, 0.125, \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0., 0., 0., 0., 0.125, \
0.166667, 0.125, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0., \
0., 0., 0., 0.125, 0.125, 0.166667, 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0.25, 0., 0., 0., 0., 0., 0., 0.166667, 0.125, 0.125, 0., \
0., 0., 0., 0., 0., 0., 0., 0.25, 0., 0., 0., 0., 0., 0., 0.125, \
0.166667, 0.125, 0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0., 0., 0., \
0., 0., 0., 0.125, 0.125, 0.166667, 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., \
0., 0., 0., 0.};

