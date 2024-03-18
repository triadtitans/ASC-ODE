#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

#include "rigid_body.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(Transformation);
PYBIND11_MAKE_OPAQUE(MassMatrix);

PYBIND11_MODULE(rigid_body, m) {
    m.doc() = "rigid body simulator";       
       
    py::class_<Transformation>(m,"Transformation")
      .def(py::init<>())
      .def("__str__", [](Transformation & t) {
        std::stringstream sstr;
        sstr << t;
        return sstr.str();
      })
      .def("setTranslation",&Transformation::setTranslation)
      .def("setRotation",&Transformation::setRotation)
      .def("asTuple",[](Transformation& t){
        return py::make_tuple(t.q_(3),t.q_(4),t.q_(5),t.q_(0),t.q_(6),t.q_(7),t.q_(8),t.q_(1),t.q_(9),t.q_(10),t.q_(11),t.q_(2));
      });

    py::class_<MassMatrix>(m,"MassMatrix")
      .def(py::init<>())
      .def("set",&MassMatrix::set);

    m.def("mass_cube",[](){return MassMatrix(MatrixView<double>(18,18,mass_matrix_data));});

    py::class_<RigidBody> (m, "RigidBody")
      .def(py::init<>())
      /*.def("__str__", [](RigidBody & rb) {
        std::stringstream sstr;
        sstr << mss;
        return sstr.str();
      })*/        
      .def_property("q", &RigidBody::getQ,&RigidBody::setQ)
      .def_property("dq", &RigidBody::getDq,&RigidBody::setDq)
      .def_property("ddq", &RigidBody::getDdq,&RigidBody::setDdq)   
      .def("setMass",&RigidBody::setMass)    
      .def("simulate",[](RigidBody& r, double tend,double steps) {r.simulate(tend,steps);}); 
    
      
    
}
