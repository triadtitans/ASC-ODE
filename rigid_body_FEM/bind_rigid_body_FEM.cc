#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

#include "rigid_body_fem.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(Transformation<>);
// PYBIND11_MAKE_OPAQUE(MassMatrix);

PYBIND11_MODULE(rigid_body_FEM, rbd) {

    // adds the BLA bindings as submodule, accessible as rigid_body.bla.Matrix etc.
    // https://github.com/pybind/pybind11/discussions/4027
    auto m = rbd.def_submodule("bla", "basic linear algebra");
    #include "bind_bla_obj.h"



    // the main bindings:
    rbd.doc() = "rigid body FEM simulator";       
       
    py::class_<Transformation<>>(rbd,"Transformation")
      .def(py::init<>())
      .def("__str__", [](Transformation<> & t) {
        std::stringstream sstr;
        sstr << t;
        return sstr.str();
      })
      .def("setTranslation",&Transformation<>::setTranslation)
      .def("setRotation",&Transformation<>::setRotation)
      .def("asTuple",[](Transformation<>& t){
        // *column-major* transformation matrix as in https://threejs.org/docs/#api/en/math/Matrix4
        // old version: return py::make_tuple(t.q_(3),t.q_(6),t.q_(9),0,t.q_(4),t.q_(7),t.q_(10),0,t.q_(5),t.q_(8),t.q_(11),0,t.q_(0),t.q_(1),t.q_(2),1);
        // new version, converts Schöberl-style ordering of Q to column-major ordering of a three.js transformation matrix:
        return py::make_tuple(t.q_(1), t.q_(5), t.q_(9), 0, t.q_(2), t.q_(6), t.q_(10), 0, t.q_(3), t.q_(7), t.q_(11), 0, t.q_(0), t.q_(4), t.q_(8), 1);
      });

    

    py::class_<RigidBody_FEM> (rbd, "RigidBody_FEM")
      .def(py::init<>())
      /*.def("__str__", [](RigidBody & rb) {
        std::stringstream sstr;
        sstr << mss;
        return sstr.str();
      })*/
      .def_property("q", &RigidBody_FEM::getQ,&RigidBody_FEM::setQ)
      .def_property("phat", &RigidBody_FEM::getPhat,&RigidBody_FEM::setPhat)
      .def_property("center",
        [](RigidBody_FEM& r){
          std::array<double,3> t;
          for (int i = 0;i<3;i++) t[i]=r.center()(i);
          return t;
        },
        [](RigidBody_FEM& r, std::array<double,3> t){
          for (int i = 0;i<3;i++) r.center()(i)=t[i];
        })
      .def_property("mass",
        [](RigidBody_FEM& r){
          return r.mass();
        },
        [](RigidBody_FEM& r, double m){
          r.mass()=m;
      })
      .def_property("inertia",
        [](RigidBody_FEM& r){
          return r.inertia();
        },
        [](RigidBody_FEM& r, Matrix<double> m){
          r.inertia()=m;
      })
      //.def("setMass", &RigidBody::setMass)
      //.def("recalcMassMatrix", &RigidBody_FEM::recalcMassMatrix)
      .def("saveState", &RigidBody_FEM::saveState)
      .def("reset", &RigidBody_FEM::reset) 
      .def("setPhat", &RigidBody_FEM::setPhat_v) 
      .def("simulate",[](RigidBody_FEM& r, double tend,double steps) {r.simulate(tend,steps);});
    //rbd.def("mass_matrix_from_inertia", &mass_matrix_from_inertia, "generates the a mass matrix from given inertia, center and mass",
            //py::arg("inertia_matrix"), py::arg("center_of_mass"), py::arg("mass"));
  
}
