#include "rigid_body_fem.h"


int main()
{
  double tend = 0.25*10;
  double steps = 5*10;
  Vector<double> q ( 12 );
  q(0)=0; q(4)=0; q(8)=0; 

  q(1)=1; q(2)=0; q(3)=0; 
  q(5)=0; q(6)=1; q(7)=0; 
  q(9)=0; q(10)=0; q(11)=1;
  
  //q(12)=0; q(13)=0;q(14)=0; q(15)=0; q(16)=0; q(17)=0; 
  Vector<double> phat(6);
  phat = 0;
  //phat(3) = 0.01;
 
  Matrix<double> inertia_matrix(3, 3);
  MatrixView<double> inertia_v (inertia_matrix);
  RigidBody_FEM rb(q,phat,1,Vec<3>{0,0,0},inertia_v);
  q(0) = 2; q(4) = 2; q(8) = 2;
  RigidBody_FEM rb2(q,phat,1,Vec<3>{0,0,0},inertia_v);
  /* RigidBody_FEM rb;
  rb.setPhat_v(3, 0.01);

  RigidBody_FEM rb2;
  rb2.setPhat_v(3, 0.01); */

  RBS_FEM rbs;
  rbs.gravity() = {0, 0, 1};
  /* rbs.bodies().push_back(rb);
  rbs.bodies().push_back(rb2); */

  Connector c1 = rbs.addBody(rb);
  Connector c3 = rbs.addBody(rb2);
  Connector c2{ConnectorType::fix, {0, 0, 3}, 0}; //= rbs.addBody(rb2);
  /* c1.pos = {1, 2, 3};
  c2.pos = {4, 5, 6}; */
  double len = Norm(rb.absolutePosOf(c1.pos) - rb2.absolutePosOf(c2.pos));

  Beam beam{0, c1, c2};
  rbs.addBeam(beam);

  simulate(rbs,tend, steps, [](int i, double t, VectorView<double> q) { 
                    std::cout<<std::fixed << "Body1 newton-iteration: " << i << " newton-error: " << std::scientific << t << std::fixed << std::endl
                      <<"\t"<< "Translation =" << q(0) << " ," << q(1) << ", "<<", " << q(2) << "} " << std::endl
                      <<"\t"<< " Rotation: " << q(3) << " ," << q(4) << ", "<<", " << q(5) << "} " << std::endl
                      <<"\t"<< "           " << q(6) << " ," << q(7) << ", "<<", " << q(8) << "} " << std::endl
                      <<"\t"<< "           " << q(9) << " ," << q(10) << ", "<<", " << q(11) << "} " << std::endl << std::endl
                      << "Body2 newton-iteration: " << i << " newton-error: " << std::scientific << t << std::fixed << std::endl
                      <<"\t"<< "Translation =" << q(0+18) << " ," << q(1+18) << ", "<<", " << q(2+18) << "} " << std::endl
                      <<"\t"<< " Rotation: " << q(3+18) << " ," << q(4+18) << ", "<<", " << q(5+18) << "} " << std::endl
                      <<"\t"<< "           " << q(6+18) << " ," << q(7+18) << ", "<<", " << q(8+18) << "} " << std::endl
                      <<"\t"<< "           " << q(9+18) << " ," << q(10+18) << ", "<<", " << q(11+18) << "} " << std::endl << std::endl; }                   
                   );
  std::cout << "}";
}
