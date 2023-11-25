#include "mass_spring.h"

int main()
{
  MassSpringSystem<2> mss;
  mss.SetGravity( {0,-9.81} );
  auto fA = mss.AddFix( { { 0.0, 0.0 } } );
  auto mA = mss.AddMass( { 1, { 1.0, 0.0 } } );
  mss.AddSpring ( { 1, 10, { fA, mA } }  );

  auto mB = mss.AddMass( { 1, { 2.0, 0.0 } } );
  mss.AddSpring ( { 1, 20, { mA, mB } } );
  
  cout << "mss: " << endl << mss << endl;

  auto mss_func = make_shared<MSS_Function<2>> (mss);



  double tend = 10;
  double dt = tend/1000;
  
  Vector<> x(2*mss.Masses().size());
  Vector<> dx(2*mss.Masses().size());  
  Vector<> ddx(2*mss.Masses().size());  

  mss.GetState (x, dx, ddx);
  
  SolveODE_Verlet(tend, dt, x, dx,  mss_func,
                  [](double t, VectorView<double> x) { cout << "t = " << t
                                                            << ", x = " << Vec<4>(x) << endl; });
}
