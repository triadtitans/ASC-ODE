#include "mass_spring.h"

int main()
{
  MassSpringSystem<2> mss;
  mss.SetGravity( {0,-9.81} );
  auto fA = mss.AddFix( { { 0.0, 0.0 } } );
  auto mA = mss.AddMass( { 1, { 1.0, 0.0 } } );
  mss.AddBeam ( { 1, { fA, mA } }  );

  auto mB = mss.AddMass( { 1, { 2.0, 0.0 } } );
  mss.AddSpring ( { 1, 20, { mA, mB } } );
  
  std::cout << "mss: " << std::endl << mss << std::endl;


  double tend = 10;
  double steps = 1000;
  
  Vector<> x(2*mss.Masses().size()+1);
  Vector<> dx(2*mss.Masses().size()+1);  
  Vector<> ddx(2*mss.Masses().size()+1);  

  auto mss_func = std::make_shared<MSS_Function<2>> (mss);
  auto mass = std::make_shared<Projector> (x.Size(),0,4);      

  mss.GetState (x, dx, ddx);
  
  SolveODE_Alpha(tend, steps,0.8, x, dx, ddx, mss_func, mass,
                   [](double t, VectorView<double> x) { std::cout << "t = " << t
                                                             << ", x = " << Vec<4>(x) << std::endl; });
}
