#include "rigid_body_fem.h"


int main()
{
  double tend = 50*2*M_PI;
  double steps = 10000;
  Vector<double> q { 18 };
  q(0)=0; q(4)=0; q(8)=0; 

  q(1)=1; q(2)=0; q(3)=0; 
  q(5)=0; q(6)=1; q(7)=0; 
  q(9)=0; q(10)=0; q(11)=1;
  
  q(12)=0; q(13)=0;q(14)=0; q(15)=0; q(16)=0; q(17)=0; 
  Vector<double> phat(6);
  phat = 0;
  phat(3) = 1;
 
  Matrix<double> inertia_matrix(3, 3);
  MatrixView<double> inertia_v (inertia_matrix);
  RigidBody rb(q,phat,1,Vec<3>{0,0,0},inertia_v);
  rb.simulate(tend,steps , [](int i, double t, VectorView<double> q) { 
                    std::cout<<std::fixed << "newton-iteration: " << i << " newton-error: " << t << std::endl
                      <<"\t"<< "Translation =" << q(0) << " ," << q(1) << ", "<<", " << q(2) << "} " << std::endl
                      <<"\t"<< " Rotation: " << q(3) << " ," << q(4) << ", "<<", " << q(5) << "} " << std::endl
                      <<"\t"<< "           " << q(6) << " ," << q(7) << ", "<<", " << q(8) << "} " << std::endl
                      <<"\t"<< "           " << q(9) << " ," << q(10) << ", "<<", " << q(11) << "} " << std::endl << std::endl; }                   
                   );
  std::cout << "}";
}
