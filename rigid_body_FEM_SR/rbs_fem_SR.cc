#include "rbs_fem_SR.h"
#include "eq_system_fem_SR.h"


int main()
{
  double tend = 0.15;
  double steps = 2;
  Vector<double> q ( 12 );
  q(0)=0; q(1)=1; q(2)=1;

  q(3)=1; q(4)=0; q(5)=0;
  q(6)=0; q(7)=1; q(8)=0;
  q(9)=0; q(10)=0; q(11)=1;

  Vector<double> q1 ( 12 );
  q1(0)=0; q1(1)=0; q1(2)=0;

  q1(3)=1; q1(4)=0; q1(5)=0;
  q1(6)=0; q1(7)=1; q1(8)=0;
  q1(9)=0; q1(10)=0; q1(11)=1;

  Vector<double> phat(6);
  phat = 0;
  //phat(3) = 0.01;

  Matrix<double> inertia_matrix = Diagonal(3, 1);
  MatrixView<double> inertia_v (inertia_matrix);
  RigidBody_FEM rb(q,phat,1,Vec<3>{0,0,0},inertia_v);
  //q(0) = 3; q(1) = 3; q(2) = 3;
  rb.recalcMassMatrix();
  RigidBody_FEM rb2(q1,phat,1,Vec<3>{0,0,0},inertia_v);
  rb2.recalcMassMatrix();
  
   
  RBS_FEM rbs;
  rbs.Gravity() = {0, 0, 9.81};

  Connector c1 = rbs.add(rb);
  Connector f1{ConnectorType::fix, {0, 0, 1}, 0}; //= rbs.addBody(rb2);
  Connector c2 = rbs.add(rb2);

  Beam bm1(c1, c2);
  rbs.add(bm1);

  Beam bm2(c1, f1);
  rbs.add(bm2);

  rbs.info_rbs();
  rbs.Energy();

  for(size_t z = 0; z < 500; z++) {
    simulate(rbs, 0.003, 1, [](int i, double t, VectorView<double> q) {
                    std::cout<<std::fixed << "Body1 newton-iteration: " << i << " newton-error: " << std::scientific << t << std::fixed << std::endl
                      <<"\t"<< "Translation =" << q(0) << " ," << q(1) << ", "<<", " << q(2) << "} " << std::endl
                      <<"\t"<< " Rotation: " << q(3) << " ," << q(4) << ", "<<", " << q(5) << "} " << std::endl
                      <<"\t"<< "           " << q(6) << " ," << q(7) << ", "<<", " << q(8) << "} " << std::endl
                      <<"\t"<< "           " << q(9) << " ," << q(10) << ", "<<", " << q(11) << "} " << std::endl << std::endl; }
                   );
  }
}
