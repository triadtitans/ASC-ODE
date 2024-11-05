#include "rbs_fem_GGL.h"
#include "eq_system_fem_GGL.h"


int main()
{
  double tend = 0.15;
  double steps = 2;
  Vector<double> q ( 12 );
  q(0)=0; q(1)=0; q(2)=1;

  q(3)=1; q(4)=0; q(5)=0;
  q(6)=0; q(7)=1; q(8)=0;
  q(9)=0; q(10)=0; q(11)=1;

  Vector<double> phat(6);
  phat = 0;

  Matrix<double> inertia_matrix = Diagonal(3, 1.0);
  MatrixView<double> inertia_v (inertia_matrix);
  RigidBody_FEM rb(q,phat,1,Vec<3>{0,0,0},inertia_v);

  rb.recalcMassMatrix();
  //RigidBody_FEM rb2(q1, phat, 1, Vec<3>{0,0,0}, inertia_v);
  //rb2.recalcMassMatrix();


  RBS_FEM rbs;
  rbs.Gravity() = {0, 9.81, 0};

  Connector c1 = rbs.add(rb);
  c1.Pos(0) = 0.5;
  Connector c2{ConnectorType::mass, Vector<double>(3), 0};
  c2.Pos(0) = -0.5;
  Connector f1{ConnectorType::fix, {0, 0, 0}, 0}; //= rbs.addBody(rb2);
  //Connector c2 = rbs.add(rb2);

  Beam bm1(f1, c1);
  rbs.add(bm1);

  //Beam bm2(f1, c2);
  //rbs.add(bm2);

  rbs.info_rbs();
  rbs.Energy();

  for( size_t i = 0; i <1 ; i++) {
    simulate(rbs,tend, steps, [](int i, double t, VectorView<double> q) {
                      std::cout<<std::fixed << "Body1 newton-iteration: " << i << " newton-error: " << std::scientific << t << std::fixed << std::endl
                        <<"\t"<< "Translation =" << q(0) << " ," << q(1) << ", "<<", " << q(2) << "} " << std::endl
                        <<"\t"<< " Rotation: " << q(3) << " ," << q(4) << ", "<<", " << q(5) << "} " << std::endl
                        <<"\t"<< "           " << q(6) << " ," << q(7) << ", "<<", " << q(8) << "} " << std::endl
                        <<"\t"<< "           " << q(9) << " ," << q(10) << ", "<<", " << q(11) << "} " << std::endl << std::endl;
                        }
                    );
    std::cout << "}";
    rbs.Energy();
  }
}
