#ifndef rbs_eq
#define rbs_eq

#include <nonlinfunc.h>
#include <math.h>
#include <matrix.h>
#include <vector.h>
#include <ode.h>
#include "../src/autodiffdiff.h"
#include "rigid_body_helper.h"
#include "rbs_fem_clean.h"
#include <cmath>
#include <chrono>


class EQRigidBody;
class EQRBS;


using namespace ASC_bla;
using namespace ASC_ode;


// system of equations for one rigid body
class EQRigidBody : public NonlinearFunction
{
  // timestep duration
  double h_;
  
  // current body transformation
  Transformation<double> Q_;
  // current momentum
  Vector<double> phatold_;
  //  current force;
  Vector<double> force_old_;
  //  current constraint
  Vector<double> vel_con_old_;

  Matrix<double> G_old_;
  // parent rigid body system
  RBS_FEM& rbs_;
  // index within bodies_ of parent rbs_
  size_t body_index_;
  
  public:
  EQRigidBody(Transformation<double> Q, Vector<double> Phat, double h, RBS_FEM& rbs, size_t body_index, Vector<double> force_old, Vector<double> vel_con_old, Matrix<double> G_old):
      h_(h), Q_(Q), phatold_(Phat), rbs_(rbs), body_index_(body_index), force_old_(force_old), vel_con_old_(vel_con_old), G_old_(G_old) {};

  size_t DimX() const  { return eq_per_body; }
  size_t DimF() const  { return eq_per_body; }

  size_t Index()  const {
    return body_index_;
  }

  void Evaluate (VectorView<double> x, VectorView<double> f) const override {
    // x (a(3), B(row_maj)(9), v_trans(3), v_skew(3), phat(6), p(6))

    // variables
    Vector<double> anew = x.Range(Index()*dim_per_body, Index()*dim_per_body + 3);
    MatrixView<double> Bnew = AsMatrix(x.Range(Index()*dim_per_body + 3, Index()*dim_per_body + 12), 3, 3);
    Vector<double> vtrans = x.Range(Index()*dim_per_body + 12, Index()*dim_per_body + 15);
    Vector<double> vskew = x.Range(Index()*dim_per_body + 15, Index()*dim_per_body + 18);
    Vector<double> phat = x.Range(Index()*dim_per_body + 18, Index()*dim_per_body + 24);
    Vector<double> p = x.Range(Index()*dim_per_body + 24, Index()*dim_per_body + 30);

    // known constants
    Vector<double> aold = Q_.getTranslation();
    Matrix<double> Bold = Q_.getRotation();
    Vector<double> force_new(dim_per_transform);
    Vector<double> vel_con_new(dim_per_state);
    Vector<double> g_con_new(dim_per_state);
    Vector<double> g_con_old(dim_per_state);
    
    rbs_.force(x, force_new, body_index_);
    rbs_.dvelocity_constraint(x,vel_con_new, body_index_);
    Matrix<double> G_new(dim_per_transform, rbs_.NumBeams());
    rbs_.G_body(x, body_index_, G_new);

    g_con_new = G_new*x.Range(dim_per_body*rbs_.NumBodies(), x.Size()).Slice(1, 2);
    g_con_old = G_old_*x.Range(dim_per_body*rbs_.NumBodies(), x.Size()).Slice(0, 2);
    //std::cout << "G_old: " << G_old_ << std::endl;
    //std::cout << "G_new: " << G_new << std::endl;
    //std::cout << "g_con_new: " << g_con_new << std::endl;
    //std::cout << "g_con_old: " << g_con_old << std::endl;

    //std::cout << "x: " << x.Range(dim_per_body * body_index_, dim_per_body * (body_index_ + 1)) << std::endl;
    //std::cout << "x2: " << x.Range(30, 42) << std::endl;
    //std::cout << "lambda 1: " << x(dim_per_body * rbs_.NumBodies()) << std::endl;
    //std::cout << "lambda 2: " << x(dim_per_body * rbs_.NumBodies() + 1) << std::endl;
    //std::cout << "x: " << x << endl;
    //std::cout << "Body: " << body_index_ << " Force old: " << force_old_ << std::endl;
    //std::cout << "Body: " << body_index_ << " Force new: " << force_new << std::endl;
    //std::cout << "Body: " << body_index_ << " velcon old: " << vel_con_old_ << std::endl;
    //std::cout << "Body: " << body_index_ << "velcon new: " << vel_con_new << std::endl;

    Matrix<double> Bhalf(3, 3);
    Bhalf = 0.5*(Bnew + Bold);

    // std::cout << vskew << "\n" << std::endl;

    // I
    f.Range(0, 3) = (1/h_)*(anew - aold) - vtrans; // - vel_con_new.Range(dim_per_transform, dim_per_transform + 3); // - dp_gv_new_.Range(0,3); // D_p(mü*D_2(C_v(q_old, phatold)))
    f.Range(3, 6) = (1/h_)*inv_hat_map(Transpose(Bhalf)*(Bnew - Bold)) - vskew; // - vel_con_new.Range(dim_per_transform + 3, dim_per_transform + 6); // dp_gv_new_.Range(3,6);

    // II
    Vector<double> vhat(6);
    vhat.Range(0, 3) = vtrans;
    vhat.Range(3, 6) = vskew;

    // Mass matrix included
    f.Range(6, 12) = rbs_.Bodies(body_index_).Mass_matrix()*vhat - p;

    //TODO: Are Bold and Bnew in the right order?
    // III - first half
    /* Vector<double> a_force(3);
    Matrix<double> B_force(3, 3); */
    
    // force<double, double> (0.5*(aold+anew), Bhalf, a_force, B_force);
    f.Range(12, 15) = (2/h_)*(p.Range(0, 3) - phatold_.Range(0, 3)) - force_old_.Range(0, 3) - g_con_old.Range(0, 3); // - vel_con_old_.Range(0, 3); // - rbs_.get_translation(dq_gv_old_, 0); // Translation(  -D_q(mü*D_1(C_v(q_old,phat_old)))  )
    f.Range(15, 18) = inv_hat_map((Transpose(Bold) * ((2/h_)*(Bhalf*hat_map(p.Range(3, 6)) - Bold*hat_map(phatold_.Range(3, 6))) - AsMatrix(force_old_.Range(3, 12), 3, 3) - AsMatrix(g_con_old.Range(3, 12), 3, 3) ) ) ); // - AsMatrix(vel_con_old_.Range(3, 12), 3, 3) ) ) ); // - rbs_.get_rotation(dq_gv_old_, 0)))); //TODO: Are Bold and Bnew in the right order?

    // III - second half
    // force<double, double> (anew, Bnew, a_force, B_force);
    f.Range(18, 21) = (2/h_)*(phat.Range(0, 3) - p.Range(0, 3)) - force_new.Range(0, 3) - g_con_new.Range(0, 3); // - vel_con_new.Range(0, 3); //- rbs_.get_translation(dq_gv_new_, 0); // Translation(  -D_q(mü*D_1(C_v(q_new,phat_old)))  )
    // std::cout << get_translation(force_new_, 0) << std::endl;
    f.Range(21, 24) = inv_hat_map((Transpose(Bnew) * ((2/h_)*(Bnew*hat_map(phat.Range(3, 6)) - Bhalf*hat_map(p.Range(3, 6))) - AsMatrix(force_new.Range(3, 12), 3 ,3) - AsMatrix(g_con_new.Range(3, 12), 3, 3) ) ) ); // - AsMatrix(vel_con_new.Range(3, 12), 3, 3) ) ) ); //- rbs_.get_rotation(dq_gv_new_, 0)))); //      TODO: Are Bold and Bnew in the right order?

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


class EQRBS : public NonlinearFunction  {
  RBS_FEM& rbs_;
  std::shared_ptr<StackedFunction_large_input> _func;
  std::vector<std::shared_ptr<EQRigidBody>> _functions; // for body equations
  size_t _numbodies_;
  size_t _num_beams;

  public:
  size_t DimX() const  { return eq_per_body*_numbodies_ + 2*_num_beams; }
  size_t DimF() const  { return  eq_per_body*_numbodies_ + 2*_num_beams; }

  EQRBS(RBS_FEM& rbs, Vector<double> state, size_t numbodies_, size_t num_beams, double h):
    rbs_(rbs), _numbodies_(numbodies_), _num_beams(num_beams)
  {
    _func = std::make_shared<StackedFunction_large_input>();

    Vector<double> x(dim_per_body*rbs.NumBodies() + 2*rbs.NumBeams());
    x = rbs.stateToX(state);
    
    // rechnen potential und force => welche parameter für 1 body =>
    for(int i=0;i<numbodies_;i++) {

      Vector<double> f(dim_per_transform);
      rbs.force(x, f, i);

      Matrix<double> G_bod(dim_per_transform, rbs_.NumBeams());
      //std::cout << "force: " << f << std::endl;
      //std::cout << "G: " << G_body << std::endl;

      Vector<double> vel_con(dim_per_state);
      rbs.dvelocity_constraint(x, vel_con, i);

      rbs_.G_body(x, i, G_bod);
      
      auto q = x.Range(i*dim_per_body,i*dim_per_body+12);
      std::cout << "q: " << q << std::endl;
      auto p = x.Range(i*dim_per_body+18,i*dim_per_body+24);
      std::cout << "p: " << p << std::endl;
      std::shared_ptr<EQRigidBody> eq = std::make_shared<EQRigidBody>(Transformation<double>(q), p, h, rbs_, i, f, vel_con, G_bod); //(Transformation<double>(q), p, h, rbs_, i);
      
      _func->addFunction(eq);
      _functions.push_back(eq);
    }
  }

  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {  
    //std::cout << "x: " << x << std::endl;
    _func->Evaluate(x,f);
    for (size_t i = 0; i < rbs_.NumBeams(); i++)  {
      f(rbs_.NumBodies()*dim_per_body + 2*i) = rbs_.g(x, i);
      f(rbs_.NumBodies()*dim_per_body + 2*i + 1) = rbs_.velocity_constraint(x, i);
    }
    
  }
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    dNumeric(*this,x,df);
  }
};



void simulate(RBS_FEM& rbs, double tend, double steps, std::function<void(int,double,VectorView<double>)> callback = nullptr ){

  Vector<double> state(dim_per_state*rbs.NumBodies() + 2*rbs.NumBeams());
  Vector<double> x(dim_per_body*rbs.NumBodies() + 2*rbs.NumBeams());

  // copy over current values 
  rbs.getState(state);
  //std::cout << state << std::endl;
  x = rbs.stateToX(state);
  // std::cout << "first x: " << x << std::endl;
  for (size_t step=0; step < steps; step++){
    std::shared_ptr<EQRBS> eq = std::make_shared<EQRBS>(rbs, state, rbs.NumBodies(), rbs.NumBeams(), tend/steps);
    // solve equation
    // Matrix<double> df(60,60);
    //Vector<double> test(30);
    //eq->Evaluate(x,test);
    //auto mv = df.Cols(0,30).Rows(0,30);
    //Matrix<double> mv(30,30);
    //eq->EvaluateDeriv(x,df);
    //std::cout << test<< std::endl;
    //std::cout << df << std::endl;
    // std::cout << "before newton" << std::endl;
    NewtonSolver(eq, x, 1e-10, 10, callback);//, [](int a,  double b, auto c){std::cout << "new_run" << std::endl;});
    //std::cout <<"l"<< x<<std::endl;
    state = rbs.xToState(x);

  
    // std::cout << state << std::endl;

    //store data into different bodies
    rbs.setState(state);
    

    /* Transformation<double> t = rbs.bodies()[1].q();
    std::cout << t << std::endl; */
  }

} 

#endif