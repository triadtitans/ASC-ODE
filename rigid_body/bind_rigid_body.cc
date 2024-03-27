#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

#include "rigid_body.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(Transformation);
// PYBIND11_MAKE_OPAQUE(MassMatrix);

PYBIND11_MODULE(rigid_body, rbd) {

    // adds the BLA bindings as submodule, accessible as rigid_body.bla.Matrix etc.
    // https://github.com/pybind/pybind11/discussions/4027
    auto m = rbd.def_submodule("bla", "basic linear algebra");
    #include "bind_bla_obj.h"



    // the main bindings:
    rbd.doc() = "rigid body simulator";       
       
    py::class_<Transformation>(rbd,"Transformation")
      .def(py::init<>())
      .def("__str__", [](Transformation & t) {
        std::stringstream sstr;
        sstr << t;
        return sstr.str();
      })
      .def("setTranslation",&Transformation::setTranslation)
      .def("setRotation",&Transformation::setRotation)
      .def("asTuple",[](Transformation& t){
        // return py::make_tuple(t.q_(3),t.q_(4),t.q_(5),t.q_(0),t.q_(6),t.q_(7),t.q_(8),t.q_(1),t.q_(9),t.q_(10),t.q_(11),t.q_(2),0,0,0,1);
        // *column-major* transformation matrix as in https://threejs.org/docs/#api/en/math/Matrix4
        // old version: return py::make_tuple(t.q_(3),t.q_(6),t.q_(9),0,t.q_(4),t.q_(7),t.q_(10),0,t.q_(5),t.q_(8),t.q_(11),0,t.q_(0),t.q_(1),t.q_(2),1);
        // new version, converts Schöberl-style ordering of Q to column-major ordering of a three.js transformation matrix:
        return py::make_tuple(t.q_(1), t.q_(5), t.q_(9), 0, t.q_(2), t.q_(6), t.q_(10), 0, t.q_(3), t.q_(7), t.q_(11), 0, t.q_(0), t.q_(4), t.q_(8), 1);
      });

    /* py::class_<MassMatrix>(rbd,"MassMatrix")
      .def(py::init<>())
      .def("set",&MassMatrix::set); */

    // rbd.def("mass_cube",[](){return MassMatrix(MatrixView<double>(18,18,mass_matrix_data));});

    // Matrix<double> MASS_CUBE (18, 18);
    // MASS_CUBE = MatrixView<double>(18, 18, mass_matrix_data);
    // rbd.attr("MASS_CUBE") = MASS_CUBE;

    py::class_<RigidBody> (rbd, "RigidBody")
      .def(py::init<>())
      /*.def("__str__", [](RigidBody & rb) {
        std::stringstream sstr;
        sstr << mss;
        return sstr.str();
      })*/        
      .def_property("q", &RigidBody::getQ,&RigidBody::setQ)
      .def_property("dq", &RigidBody::getDq,&RigidBody::setDq)
      .def_property("ddq", &RigidBody::getDdq,&RigidBody::setDdq)   
      .def("setMass", &RigidBody::setMass)
      .def("saveState", &RigidBody::saveState)
      .def("reset", &RigidBody::reset)  
      .def("simulate",[](RigidBody& r, double tend,double steps) {r.simulate(tend,steps);});
    

    rbd.def("mass_matrix_from_inertia", &mass_matrix_from_inertia, "generates the a mass matrix from given inertia, center and mass",
            py::arg("inertia_matrix"), py::arg("center_of_mass"), py::arg("mass"));
    
}
