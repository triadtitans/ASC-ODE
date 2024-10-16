#ifndef rbs_eq
#define rbs_eq

#include <nonlinfunc.h>
#include <math.h>
#include <matrix.h>
#include <vector.h>
#include <ode.h>
#include "../src/autodiffdiff.h"
#include "rigid_body_helper_GGL.h"
#include "rbs_fem_GGL.h"
#include <cmath>
#include <chrono>


class EQRigidBody;
class EQRBS;


using namespace ASC_bla;
using namespace ASC_ode;


// system of equations for one rigid body
class EQRigidBody : public NonlinearFunction
{
  //  timestep duration
  double h_;
  //  current body transformation
  Transformation<double> Q_;
  //  current momentum
  Vector<double> phatold_;
  //  current force;
  Vector<double> force_old_;
  //  current constraint
  Vector<double> vel_con_old_;
  //  parent rigid body system
  RBS_FEM& rbs_;
  //  index within bodies_ of parent rbs_
  size_t body_index_;
  
  public:
  EQRigidBody(Transformation<double> Q, Vector<double> Phat, double h, RBS_FEM& rbs, size_t body_index, Vector<double> force_old, Vector<double> vel_con_old):
      h_(h), Q_(Q), phatold_(Phat), rbs_(rbs), body_index_(body_index), force_old_(force_old), vel_con_old_(vel_con_old) {};

  size_t DimX() const  { return eq_per_body; }
  size_t DimF() const  { return eq_per_body; }

  //  get-method body index
  size_t Index()  const {
    return body_index_;
  }

  void Evaluate (VectorView<double> x, VectorView<double> f) const override {
    //  x (a(3), B(row_maj)(9), v_trans(3), v_skew(3), phat(6), p(6))

    //  extract variables
    Vector<double> anew = x.Range(Index()*dim_per_body, Index()*dim_per_body + 3);
    MatrixView<double> Bnew = AsMatrix(x.Range(Index()*dim_per_body + 3, Index()*dim_per_body + 12), 3, 3);
    Vector<double> vtrans = x.Range(Index()*dim_per_body + 12, Index()*dim_per_body + 15);
    Vector<double> vskew = x.Range(Index()*dim_per_body + 15, Index()*dim_per_body + 18);
    Vector<double> phat = x.Range(Index()*dim_per_body + 18, Index()*dim_per_body + 24);
    Vector<double> p = x.Range(Index()*dim_per_body + 24, Index()*dim_per_body + 30);

    //  get known constants
    Vector<double> aold = Q_.getTranslation();
    Matrix<double> Bold = Q_.getRotation();

    //  prepare velocity constraint and force
    Vector<double> force_new(dim_per_transform);
    Vector<double> vel_con_new(dim_per_state);

    //  calculate velocity constraint and force
    rbs_.force(x, force_new, body_index_);
    rbs_.dvelocity_constraint(x,vel_con_new, body_index_);

    Matrix<double> Bhalf(3, 3);
    Bhalf = 0.5*(Bnew + Bold);

    //  I
    f.Range(0, 3) = (1/h_)*(anew - aold) - vtrans - vel_con_new.Range(dim_per_transform, dim_per_transform + 3);
    f.Range(3, 6) = (1/h_)*inv_hat_map(Transpose(Bhalf)*(Bnew - Bold)) - vskew - vel_con_new.Range(dim_per_transform + 3, dim_per_transform + 6);

    //  II
    Vector<double> vhat(6);
    vhat.Range(0, 3) = vtrans;
    vhat.Range(3, 6) = vskew;

    //  Mass matrix included
    f.Range(6, 12) = rbs_.Bodies(body_index_).Mass_matrix()*vhat - p;

    //  III - first half
    f.Range(12, 15) = (2/h_)*(p.Range(0, 3) - phatold_.Range(0, 3)) - force_old_.Range(0, 3) - vel_con_old_.Range(0, 3);
    f.Range(15, 18) = inv_hat_map((Transpose(Bold) * ((2/h_)*(Bhalf*hat_map(p.Range(3, 6)) - Bold*hat_map(phatold_.Range(3, 6))) - AsMatrix(force_old_.Range(3, 12), 3, 3) - AsMatrix(vel_con_old_.Range(3, 12), 3, 3) ) ) );

    //  III - second half
    f.Range(18, 21) = (2/h_)*(phat.Range(0, 3) - p.Range(0, 3)) - force_new.Range(0, 3) - vel_con_new.Range(0, 3);
    f.Range(21, 24) = inv_hat_map((Transpose(Bnew) * ((2/h_)*(Bnew*hat_map(phat.Range(3, 6)) - Bhalf*hat_map(p.Range(3, 6))) - AsMatrix(force_new.Range(3, 12), 3 ,3) - AsMatrix(vel_con_new.Range(3, 12), 3, 3) ) ) );

    //  Bnew must be orthonormal
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

    df = 0;

    // setup vector to take derivative
    Vector<AutoDiffDiff<dim_per_body, double>> x_diff = x;
    Vector<AutoDiffDiff<dim_per_body, double>> f(dim_per_body);

    size_t last_bd_index = Index();
    size_t curr_bd_index = Index();
    size_t count = 0;

    do {
      curr_bd_index = rbs_.Bodies(Index()).Connected_Bodies(count);
      
      // derivative after current body
      for(size_t i = 0; i < dim_per_body; i++) {
        x_diff(last_bd_index * dim_per_body + i).DValue(i) = 0;
        x_diff(curr_bd_index * dim_per_body + i).DValue(i) = 1;
      }

      // extract variables
      Vector<AutoDiffDiff<dim_per_body, double>> anew = x_diff.Range(Index()*dim_per_body, Index()*dim_per_body + 3);
      Matrix<AutoDiffDiff<dim_per_body, double>> Bnew = AsMatrix(x_diff.Range(Index()*dim_per_body + 3, Index()*dim_per_body + 12), 3, 3);
      Vector<AutoDiffDiff<dim_per_body, double>> vtrans = x_diff.Range(Index()*dim_per_body + 12, Index()*dim_per_body + 15);
      Vector<AutoDiffDiff<dim_per_body, double>> vskew = x_diff.Range(Index()*dim_per_body + 15, Index()*dim_per_body + 18);
      Vector<AutoDiffDiff<dim_per_body, double>> phat = x_diff.Range(Index()*dim_per_body + 18, Index()*dim_per_body + 24);
      Vector<AutoDiffDiff<dim_per_body, double>> p = x_diff.Range(Index()*dim_per_body + 24, Index()*dim_per_body + 30);

      //  get known constants
      Vector<double> aold = Q_.getTranslation();
      Matrix<double> Bold = Q_.getRotation();

      //  prepare velocity constraint and force
      Vector<AutoDiffDiff<dim_per_body, double>> force_new(dim_per_transform);
      Vector<AutoDiffDiff<dim_per_body, double>> vel_con_new(dim_per_state);

      //  calculate velocity constraint and force
      rbs_.force(x_diff, force_new, body_index_);
      rbs_.dvelocity_constraint(x_diff, vel_con_new, body_index_);

      Matrix<AutoDiffDiff<dim_per_body, double>> Bhalf(3, 3);
      Bhalf = 0.5*(Bnew + Bold);

      // I
      f.Range(0, 3) = (1/h_)*(anew - aold) - vtrans - vel_con_new.Range(dim_per_transform, dim_per_transform + 3);
      f.Range(3, 6) = (1/h_)*inv_hat_map(Transpose(Bhalf)*(Bnew - Bold)) - vskew - vel_con_new.Range(dim_per_transform + 3, dim_per_transform + 6);

      // II
      Vector<AutoDiffDiff<dim_per_body, double>> vhat(6);
      vhat.Range(0, 3) = vtrans;
      vhat.Range(3, 6) = vskew;

      // Mass matrix included
      f.Range(6, 12) = rbs_.Bodies(body_index_).Mass_matrix()*vhat - p;
      
      // III - first half
      f.Range(12, 15) = (2/h_)*(p.Range(0, 3) - phatold_.Range(0, 3)) - force_old_.Range(0, 3) - vel_con_old_.Range(0, 3);
      f.Range(15, 18) = inv_hat_map((Transpose(Bold) * ((2/h_)*(Bhalf*hat_map(p.Range(3, 6)) - Bold*hat_map(phatold_.Range(3, 6))) - AsMatrix(force_old_.Range(3, 12), 3, 3) - AsMatrix(vel_con_old_.Range(3, 12), 3, 3) ) ) );

      // III - second half
      f.Range(18, 21) = (2/h_)*(phat.Range(0, 3) - p.Range(0, 3)) - force_new.Range(0, 3) - vel_con_new.Range(0, 3);
      f.Range(21, 24) = inv_hat_map((Transpose(Bnew) * ((2/h_)*(Bnew*hat_map(phat.Range(3, 6)) - Bhalf*hat_map(p.Range(3, 6))) - AsMatrix(force_new.Range(3, 12), 3 ,3) - AsMatrix(vel_con_new.Range(3, 12), 3, 3) ) ) );

      // Bnew must be orthonormal
      Matrix<double> eye = Diagonal(3, 1);
      auto c = Transpose(Bnew)* Bnew + (-1)*eye;
      f(0 + 24)=c(0, 0);
      f(1 + 24)=c(1, 0);
      f(2 + 24)=c(2, 0);
      f(3 + 24)=c(1, 1);
      f(4 + 24)=c(2, 1);
      f(5 + 24)=c(2, 2);


      count += 1;
      last_bd_index = curr_bd_index;

      MatrixView<double> view = df.Cols(curr_bd_index* dim_per_body, dim_per_body);

      for (size_t i = 0; i < dim_per_body; i++) {
        view.Row(i) = f(i).DValue_vec();
      }

    } while (count < rbs_.Bodies(Index()).Connected_Bodies().size());


    Vector<AutoDiffDiff<2, double>> xb_diff = x;
    Vector<AutoDiffDiff<2, double>> fb(dim_per_body);

    size_t curr_beam_index = 0;
    size_t prev_beam_index = 0;
    count = 0;
    
    while (count < rbs_.Bodies(Index()).Beams().size()) {
      fb = 0;
      curr_beam_index = rbs_.Bodies(Index()).Beams(count);
      
      //  derivative after current body
      xb_diff(rbs_.NumBeams() * dim_per_body + 2 * prev_beam_index).DValue(0) = 0;
      xb_diff(rbs_.NumBeams() * dim_per_body + 2 * curr_beam_index).DValue(0) = 1;

      xb_diff(rbs_.NumBeams() * dim_per_body + 2 * prev_beam_index + 1).DValue(1) = 0;
      xb_diff(rbs_.NumBeams() * dim_per_body + 2 * curr_beam_index + 1).DValue(1) = 1;

      //  extract variables
      Vector<AutoDiffDiff<2, double>> anew = xb_diff.Range(Index()*dim_per_body, Index()*dim_per_body + 3);
      Matrix<AutoDiffDiff<2, double>> Bnew = AsMatrix(xb_diff.Range(Index()*dim_per_body + 3, Index()*dim_per_body + 12), 3, 3);
      Vector<AutoDiffDiff<2, double>> vtrans = xb_diff.Range(Index()*dim_per_body + 12, Index()*dim_per_body + 15);
      Vector<AutoDiffDiff<2, double>> vskew = xb_diff.Range(Index()*dim_per_body + 15, Index()*dim_per_body + 18);
      Vector<AutoDiffDiff<2, double>> phat = xb_diff.Range(Index()*dim_per_body + 18, Index()*dim_per_body + 24);
      Vector<AutoDiffDiff<2, double>> p = xb_diff.Range(Index()*dim_per_body + 24, Index()*dim_per_body + 30);

      //  get known constants
      Vector<double> aold = Q_.getTranslation();
      Matrix<double> Bold = Q_.getRotation();

      //  prepare velocity constraint and force
      Vector<AutoDiffDiff<2, double>> force_new(dim_per_transform);
      Vector<AutoDiffDiff<2, double>> vel_con_new(dim_per_state);
      
      //  calculate prepare velocity constraint and force
      rbs_.force(xb_diff, force_new, body_index_);
      rbs_.dvelocity_constraint(xb_diff, vel_con_new, body_index_);

      Matrix<AutoDiffDiff<2, double>> Bhalf(3, 3);
      Bhalf = 0.5*(Bnew + Bold);

      //  I
      fb.Range(0, 3) = (1/h_)*(anew - aold) - vtrans - vel_con_new.Range(dim_per_transform, dim_per_transform + 3);
      fb.Range(3, 6) = (1/h_)*inv_hat_map(Transpose(Bhalf)*(Bnew - Bold)) - vskew - vel_con_new.Range(dim_per_transform + 3, dim_per_transform + 6);
      
      //  III - first half
      fb.Range(12, 15) = (2/h_)*(p.Range(0, 3) - phatold_.Range(0, 3)) - force_old_.Range(0, 3) - vel_con_old_.Range(0, 3);
      fb.Range(15, 18) = inv_hat_map((Transpose(Bold) * ((2/h_)*(Bhalf*hat_map(p.Range(3, 6)) - Bold*hat_map(phatold_.Range(3, 6))) - AsMatrix(force_old_.Range(3, 12), 3, 3) - AsMatrix(vel_con_old_.Range(3, 12), 3, 3) ) ) );

      //  III - second half
      fb.Range(18, 21) = (2/h_)*(phat.Range(0, 3) - p.Range(0, 3)) - force_new.Range(0, 3) - vel_con_new.Range(0, 3);
      fb.Range(21, 24) = inv_hat_map((Transpose(Bnew) * ((2/h_)*(Bnew*hat_map(phat.Range(3, 6)) - Bhalf*hat_map(p.Range(3, 6))) - AsMatrix(force_new.Range(3, 12), 3 ,3) - AsMatrix(vel_con_new.Range(3, 12), 3, 3) ) ) );

      count += 1;
      prev_beam_index = curr_beam_index;

      MatrixView<double> view = df.Cols(rbs_.NumBodies() * dim_per_body + 2 * curr_beam_index, 2);

      for (size_t i = 0; i < dim_per_body; i++) {
        view.Row(i) = fb(i).DValue_vec();
      }
    }    
    
  }
};


class EQRBS : public NonlinearFunction  {

  //  system
  RBS_FEM& rbs_;

  //  stacked function of single body functions
  std::shared_ptr<StackedFunction_large_input> func_;
  //  single body functions
  std::vector<std::shared_ptr<EQRigidBody>> functions_; // for body equations
  //  number of bodies in the system
  size_t _numbodies_;
  //  number of beams in the system
  size_t _num_beams;

  public:
  size_t DimX() const  { return eq_per_body*_numbodies_ + 2*_num_beams; }
  size_t DimF() const  { return  eq_per_body*_numbodies_ + 2*_num_beams; }

  EQRBS(RBS_FEM& rbs, Vector<double> state, size_t numbodies_, size_t num_beams, double h):
    rbs_(rbs), _numbodies_(numbodies_), _num_beams(num_beams)
  {

    //  special stacked function in whick each function takes the whole input vector
    func_ = std::make_shared<StackedFunction_large_input>();

    //  define full size vector with state and all temporary variables needed for a timestep
    Vector<double> x(dim_per_body*rbs.NumBodies() + 2*rbs.NumBeams());

    //  initialize x
    x = rbs.stateToX(state);
    
    for(size_t i=0;i<numbodies_;i++) {
      
      //  calculate force at start of time step (not done in the body equation for performance)
      Vector<double> f(dim_per_transform);
      rbs.force(x, f, i);

      //  calculate impact of beam constraints for equations at start of time step (not done in the body equation for performance)
      Vector<double> vel_con(dim_per_state);
      rbs.dvelocity_constraint(x, vel_con, i);
      
      //  get state variables
      auto q = x.Range(i*dim_per_body,i*dim_per_body+12);
      auto p = x.Range(i*dim_per_body+18,i*dim_per_body+24);

      //  set up Equation for single body
      std::shared_ptr<EQRigidBody> eq = std::make_shared<EQRigidBody>(Transformation<double>(q), p, h, rbs_, i, f, vel_con); //(Transformation<double>(q), p, h, rbs_, i);
      
      //  add function to stacked function
      func_->addFunction(eq);
      functions_.push_back(eq);
    }
  }

  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {  
    //  |   eq body one   |
    //  |   eq body two   |
    //  |       .
    //  |       .
    //  |   eq body n     |
    //  |  constr beam 1  |
    //  |  constr beam 2  |
    //  |        .        |
    //  |        .        |
    //  |  constr beam k  |

    //  Evaluate all body equations
    func_->Evaluate(x,f);
    //std::cout << "here!" << std::endl;
    
    //  add the constraints at the end
    for (size_t i = 0; i < rbs_.NumBeams(); i++)  {
      f(rbs_.NumBodies()*dim_per_body + 2*i) = rbs_.g(x, i);
      f(rbs_.NumBodies()*dim_per_body + 2*i + 1) = rbs_.velocity_constraint(x, i);
    }
    //std::cout << "again here!" << std::endl;
    
    
  }
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {

    //  extract the all rows with body equations
    MatrixView<double> current_block = df.Rows(0, dim_per_body*rbs_.NumBodies());

    //  derivative of all body equations
    func_->EvaluateDeriv(x, current_block);


    //  set up derivative of constraint equations
    Vector<AutoDiffDiff<2*dim_per_state, double>> x_diff = x;
    
    size_t prev_bm_index = 0;
    
    //  take derivative of each beam constraints
    for (size_t j = 0; j < rbs_.NumBeams(); j++)  {

      Beam curr_bm = rbs_.Beams(j);
      Beam prev_bm = rbs_.Beams(prev_bm_index); 
      
      //  each constraint only depends on two bodys - setup of autodiff vector with values for those two bodies
      for (size_t i= 0; i < dim_per_transform; i++) {
        x_diff(prev_bm.Body_index_a()*dim_per_body + i).DValue(i) = 0;
        x_diff(prev_bm.Body_index_b()*dim_per_body + i).DValue(dim_per_state + i) = 0;
        x_diff(curr_bm.Body_index_a()*dim_per_body + i).DValue(i) = 1;
        x_diff(curr_bm.Body_index_b()*dim_per_body + i).DValue(dim_per_state + i) = 1;
      }

      for (size_t i = 0; i < dim_per_state - dim_per_transform; i++)  {
        x_diff(prev_bm.Body_index_a()*dim_per_body + 18 + i).DValue(dim_per_transform + i) = 0;
        x_diff(prev_bm.Body_index_b()*dim_per_body + 18 + i).DValue(dim_per_state + dim_per_transform + i) = 0;
        x_diff(curr_bm.Body_index_a()*dim_per_body + 18 + i).DValue(dim_per_transform + i) = 1;
        x_diff(curr_bm.Body_index_b()*dim_per_body + 18 + i).DValue(dim_per_state + dim_per_transform + i) = 1;
      }

      //  evaluate both equations
      Vector<double> res_g = rbs_.g(x_diff, curr_bm.Index()).DValue_vec();
      Vector<double> res_vel_con = rbs_.velocity_constraint(x_diff, curr_bm.Index()).DValue_vec();

      //  fill result into jacobi matrix
      for (size_t i= 0; i < dim_per_transform; i++) {
        df(rbs_.NumBeams() * dim_per_body + 2 * j, curr_bm.Body_index_a()*dim_per_body + i) = res_g(i);
        df(rbs_.NumBeams() * dim_per_body + 2 * j, curr_bm.Body_index_b()*dim_per_body + i) = res_g(dim_per_state + i);

        df(rbs_.NumBeams() * dim_per_body + 2 * j + 1, curr_bm.Body_index_a()*dim_per_body + i) = res_vel_con(i);
        df(rbs_.NumBeams() * dim_per_body + 2 * j + 1, curr_bm.Body_index_b()*dim_per_body + i) = res_vel_con(dim_per_state + i);
      }

      for (size_t i = 0; i < dim_per_state - dim_per_transform; i++)  {
        df(rbs_.NumBeams() * dim_per_body + 2 * j, curr_bm.Body_index_a()*dim_per_body + 18 + i) = res_g(dim_per_transform + i);
        df(rbs_.NumBeams() * dim_per_body + 2 * j, curr_bm.Body_index_b()*dim_per_body + 18 + i) = res_g(dim_per_state + dim_per_transform + i);

        df(rbs_.NumBeams() * dim_per_body + 2 * j + 1, curr_bm.Body_index_a()*dim_per_body + 18 + i) = res_vel_con(dim_per_transform + i);
        df(rbs_.NumBeams() * dim_per_body + 2 * j + 1, curr_bm.Body_index_b()*dim_per_body + 18 + i) = res_vel_con(dim_per_state + dim_per_transform + i);
      }
    }
    
  }
};



void simulate(RBS_FEM& rbs, double tend, double steps, std::function<void(int,double,VectorView<double>)> callback = nullptr ){

  Vector<double> state(dim_per_state*rbs.NumBodies() + 2*rbs.NumBeams());
  Vector<double> q(dim_per_body*rbs.NumBodies() + 2*rbs.NumBeams());

  // get current state of system
  rbs.getState(state);
  
  // initialize vector with all temporary variables
  q = rbs.stateToX(state);
  //std::cout << "x: " << x << std::endl;
  //std::cout << "state: " << state << std::endl;
  
  for (size_t step=0; step < steps; step++){
    std::shared_ptr<EQRBS> eq = std::make_shared<EQRBS>(rbs, state, rbs.NumBodies(), rbs.NumBeams(), tend/steps);
   
    NewtonSolver(eq, q, 1e-10, 10, callback);
    
    // extract only state for propagation
    state = rbs.xToState(q);

    //store data into different bodies
    rbs.setState(state);
  }
  //std::cout << "finished" << std::endl;
  /*
  std::cout<<std::fixed << "Body1 newton-iteration: " << " newton-error: " << std::scientific << std::fixed << std::endl
                      <<"\t"<< "Translation =" << q(0) << " ," << q(1) << ", "<<", " << q(2) << "} " << std::endl
                      <<"\t"<< " Rotation: " << q(3) << " ," << q(4) << ", "<<", " << q(5) << "} " << std::endl
                      <<"\t"<< "           " << q(6) << " ," << q(7) << ", "<<", " << q(8) << "} " << std::endl
                      <<"\t"<< "           " << q(9) << " ," << q(10) << ", "<<", " << q(11) << "} " << std::endl << std::endl
                      << "Body2 newton-iteration: " << " newton-error: " << std::scientific << std::fixed << std::endl
                      <<"\t"<< "Translation =" << q(0+30) << " ," << q(1+30) << ", "<<", " << q(2+30) << "} " << std::endl
                      <<"\t"<< " Rotation: " << q(3+30) << " ," << q(4+30) << ", "<<", " << q(5+30) << "} " << std::endl
                      <<"\t"<< "           " << q(6+30) << " ," << q(7+30) << ", "<<", " << q(8+30) << "} " << std::endl
                      <<"\t"<< "           " << q(9+30) << " ," << q(10+30) << ", "<<", " << q(11+30) << "} " << std::endl << std::endl;
  rbs.info_rbs();
  */
} 

#endif