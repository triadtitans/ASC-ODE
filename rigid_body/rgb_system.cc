#include "rgb_system.h"

int main(){
  double tend = 50*2*M_PI;
  double steps = 10000;
  Vector<double> q { 18 };
  q(0)=0; q(4)=0; q(8)=0; 

  q(1)=1; q(2)=0; q(3)=0; 
  q(5)=0; q(6)=1; q(7)=0; 
  q(9)=0; q(10)=0; q(11)=1; 
  
  q(12)=0; q(13)=0;q(14)=0; q(15)=0; q(16)=0; q(17)=0; 
  Vector<double> dq { 18 };
  dq(0)=0.00001;
  dq(10)=0.001;
  dq(7)=-0.001;
 
  Vector<double> ddq { 18 };
  MatrixView<double> mass_matrix(18,18,mass_matrix_data);
  RigidBody rb1(mass_matrix,q,dq,ddq);

  q(0)=0; q(4)=0; q(8)=0; 

  q(1)=1; q(2)=0; q(3)=0; 
  q(5)=0; q(6)=1; q(7)=0; 
  q(9)=0; q(10)=0; q(11)=1; 
  
  q(12)=0; q(13)=0;q(14)=0; q(15)=0; q(16)=0; q(17)=0; 
  Vector<double> dq { 18 };
  dq(0)=0.00001;
  dq(10)=0.001;
  dq(7)=-0.001;
 
  Vector<double> ddq { 18 };
  MatrixView<double> mass_matrix(18,18,mass_matrix_data);
  RigidBody rb2(mass_matrix,q,dq,ddq);

  RBSystem sys;
  sys.addBody(rb1);
  sys.addBody(rb2);

  sys.simulate(tend,steps, [](double t, VectorView<double> q) { 
                    std::cout<<std::fixed << t << ": Translation =" << q(0) << " ," << q(4) << ", "<<", " << q(8) << "} " << std::endl
                      <<"\t"<< " Rotation: " << q(1) << " ," << q(2) << ", "<<", " << q(3) << "} " << std::endl
                      <<"\t"<< "           " << q(5) << " ," << q(6) << ", "<<", " << q(7) << "} " << std::endl
                      <<"\t"<< "           " << q(9) << " ," << q(10) << ", "<<", " << q(11) << "} " << std::endl
                      << ": Translation =" << q.Range(18,36)(0) << " ," << q.Range(18,36)(4) << ", "<<", " << q.Range(18,36)(8) << "} " << std::endl
                      <<"\t"<< " Rotation: " << q.Range(18,36)(1) << " ," << q.Range(18,36)(2) << ", "<<", " << q.Range(18,36)(3) << "} " << std::endl
                      <<"\t"<< "           " << q.Range(18,36)(5) << " ," << q.Range(18,36)(6) << ", "<<", " << q.Range(18,36)(7) << "} " << std::endl
                      <<"\t"<< "           " << q.Range(18,36)(9) << " ," << q.Range(18,36)(10) << ", "<<", " << q.Range(18,36)(11) << "} " << std::endl; }    
}