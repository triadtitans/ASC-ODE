#include "rbs_fem_GGL.h"
#include "eq_system_fem_GGL.h"


int main()
{
  double tend = 0.15/2;
  double steps = 500;
  Vector<double> q ( 12 );
  q(0)=0; q(1)=1; q(2)=1;

  q(3)=1; q(4)=0; q(5)=0;
  q(6)=0; q(7)=1; q(8)=0;
  q(9)=0; q(10)=0; q(11)=1;

  //q(12)=0; q(13)=0;q(14)=0; q(15)=0; q(16)=0; q(17)=0;
  Vector<double> phat(6);
  phat = 0;
  //phat(3) = 0.01;

  Matrix<double> inertia_matrix = Diagonal(3, 1.0);
  MatrixView<double> inertia_v (inertia_matrix);
  RigidBody_FEM rb(q,phat,1,Vec<3>{0,0,0},inertia_v);
  rb.recalcMassMatrix();

  RBS_FEM rbs;
  rbs.Gravity() = {0, -9.81, 0};

  Connector c1 = rbs.add(rb);
  Connector c2{ConnectorType::fix, {0, 0, 2}, 0}; //= rbs.addBody(rb2);

  Beam bm1(c1, c2);
  rbs.add(bm1);


  rbs.info_rbs();

  //Spring spring(len, 0.1, c1, c2);
  //rbs.add(spring);

  //Spring spring2(len2, 0.1, c1, c3);
  //rbs.add(spring2);
  /*
  for (size_t i: rbs.Bodies()[1].Springs()) {
    std::cout << i << endl;
  }
  for (size_t i: rbs.Bodies()[0].Springs()) {
    std::cout << i << endl;
  }
  */
  rbs.Energy();

  for( size_t i = 0; i <100; i++) {
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
