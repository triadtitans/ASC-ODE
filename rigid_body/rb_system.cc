#include "rb_system.h"

int main(){


  double tend = 3;
  double steps = 1000;
  Vector<double> q { 18 };
  q(0)=0.0000; q(4)=0; q(8)=0; 

  q(1)=1; q(2)=0; q(3)=0; 
  q(5)=0; q(6)=1; q(7)=0; 
  q(9)=0; q(10)=0; q(11)=1; 
  
  q(12)=0; q(13)=0;q(14)=0; q(15)=0; q(16)=0; q(17)=0; 
  Vector<double> dq { 18 };
  dq(0)=0.1;
  dq(10)=0.000;
  dq(7)=-0.000;
 
  Vector<double> ddq { 18 };
  MatrixView<double> inertia_matrix(3,3,inertia_matrix_data);
  RigidBody rb1(q,dq,ddq,1,Vec<3>{0,0,0},inertia_matrix);
  q=0;dq=0;

  q(0)=3; q(4)=3; q(8)=3; 

  q(1)=1; q(2)=0; q(3)=0; 
  q(5)=0; q(6)=1; q(7)=0; 
  q(9)=0; q(10)=0; q(11)=1; 
  
  q(12)=0; q(13)=0;q(14)=0; q(15)=0; q(16)=0; q(17)=0; 
  dq(0)=0.05;
  dq(3)=0.000;
  dq(9)=-0.000;
 
  RigidBody rb2(q,dq,ddq,1,Vec<3>{0,0,0},inertia_matrix);

  RBSystem sys;
  sys.gravity() = Vec<3>{0,0,9.8};
  auto c1 = sys.addBody(rb1);
  auto c2 = sys.addBody(rb2);
  c1.pos = {1,1,1};
  c2.pos = {-1,-1,-1};
  auto fix = sys.addFix();
  auto t = sys.connectorPos(c1);
  //Beam b = {Norm(sys.connectorPos(c1)-sys.connectorPos(c2)),c1,c2};
  //sys.addBeam(b);

  Spring s = {1.7320508075688772,0.1,c1,c2};
  sys.addSpring(s);
  sys.simulate(tend,steps, [](double t, VectorView<double> q) { 
                    std::cout<<std::fixed << t << ": Translation =" << q(0) << " ," << q(4) << ", "<<", " << q(8) << "} " << std::endl
                      <<"\t"<< " Rotation: " << q(1) << " ," << q(2) << ", "<<", " << q(3) << "} " << std::endl
                      <<"\t"<< "           " << q(5) << " ," << q(6) << ", "<<", " << q(7) << "} " << std::endl
                      <<"\t"<< "           " << q(9) << " ," << q(10) << ", "<<", " << q(11) << "} " << std::endl
                      << ": Translation =" << q.Range(18,36)(0) << " ," << q.Range(18,36)(4) << ", "<<", " << q.Range(18,36)(8) << "} " << std::endl
                      <<"\t"<< " Rotation: " << q.Range(18,36)(1) << " ," << q.Range(18,36)(2) << ", "<<", " << q.Range(18,36)(3) << "} " << std::endl
                      <<"\t"<< "           " << q.Range(18,36)(5) << " ," << q.Range(18,36)(6) << ", "<<", " << q.Range(18,36)(7) << "} " << std::endl
                      <<"\t"<< "           " << q.Range(18,36)(9) << " ," << q.Range(18,36)(10) << ", "<<", " << q.Range(18,36)(11) << "} " << std::endl; });    
}