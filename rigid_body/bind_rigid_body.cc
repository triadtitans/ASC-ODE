#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

#include "rb_system.h"

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
        // *column-major* transformation matrix as in https://threejs.org/docs/#api/en/math/Matrix4
        // old version: return py::make_tuple(t.q_(3),t.q_(6),t.q_(9),0,t.q_(4),t.q_(7),t.q_(10),0,t.q_(5),t.q_(8),t.q_(11),0,t.q_(0),t.q_(1),t.q_(2),1);
        // new version, converts Schöberl-style ordering of Q to column-major ordering of a three.js transformation matrix:
        return py::make_tuple(t.q_(1), t.q_(5), t.q_(9), 0, t.q_(2), t.q_(6), t.q_(10), 0, t.q_(3), t.q_(7), t.q_(11), 0, t.q_(0), t.q_(4), t.q_(8), 1);
      });

    /* py::class_<MassMatrix>(rbd,"MassMatrix")
      .def(py::init<>())
      .def("set",&MassMatrix::set); */

    // Matrix<double> MASS_CUBE (18, 18);
    // MASS_CUBE = MatrixView<double>(18, 18, mass_matrix_data);
    // rbd.attr("MASS_CUBE") = MASS_CUBE;

    py::class_<Connector>(rbd,"Connector")
      .def_property("pos",[](Connector& c){return py::make_tuple(c.pos(0),c.pos(1),c.pos(2));},
                          [](Connector& c, std::array<double,3> t){c.pos(0)=t[0];c.pos(1)=t[1];c.pos(2)=t[2];})
      .def_property_readonly("type",[](Connector& c){return c.t == ConnectorType::mass ? 0 : 1 ;});

    py::class_<Beam>(rbd,"Beam")
      .def(py::init<>([](Connector a, Connector b, double length){return Beam{length,a,b};}))
      .def_property_readonly("length", [](Beam& b){return b.length;})
      .def_property_readonly("connectorA", [](Beam& b){return b.a;})
      .def_property_readonly("connectorB",[](Beam& b){return b.b;});

    py::class_<Spring>(rbd,"Spring")
      .def(py::init<>([](Connector a, Connector b, double length, double stiffness){return Spring{length,stiffness,a,b};}))
      .def_property_readonly("length", [](Spring& b){return b.length;})
      .def_property_readonly("stiffness", [](Spring& b){return b.length;})
      .def_property_readonly("connectorA", [](Spring& b){return b.a;})
      .def_property_readonly("connectorB",[](Spring& b){return b.b;});

    

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
      .def_property("center",
        [](RigidBody& r){
          std::array<double,3> t;
          for (int i = 0;i<3;i++) t[i]=r.center()(i);
          return t;
        },
        [](RigidBody& r, std::array<double,3> t){
          for (int i = 0;i<3;i++) r.center()(i)=t[i];
        })
      .def_property("mass",
        [](RigidBody& r){
          return r.mass();
        },
        [](RigidBody& r, double m){
          r.mass()=m;
      })
      .def_property("inertia",
        [](RigidBody& r){
          return r.inertia();
        },
        [](RigidBody& r, Matrix<double> m){
          r.inertia()=m;
      })
      //.def("setMass", &RigidBody::setMass)
      .def("recalcMassMatrix", &RigidBody::recalcMassMatrix)
      .def("saveState", &RigidBody::saveState)
      .def("reset", &RigidBody::reset)  
      .def("simulate",[](RigidBody& r, double tend,double steps) {r.simulate(tend,steps);});
    rbd.def("mass_matrix_from_inertia", &mass_matrix_from_inertia, "generates the a mass matrix from given inertia, center and mass",
            py::arg("inertia_matrix"), py::arg("center_of_mass"), py::arg("mass"));


    py::class_<RBSystem> (rbd, "RBSystem")
      .def(py::init<>())
      .def("addBody",&RBSystem::addBody)
      .def("addBeam",&RBSystem::addBeam)
      .def("addSpring",&RBSystem::addSpring)
      .def("addFix",&RBSystem::addFix)
      .def("simulate", [](RBSystem& sys,double tend, double steps){sys.simulate(tend,steps);})
      .def("bodies", &RBSystem::bodies)
      .def_property("gravity",
        [](RBSystem& r){
          std::array<double,3> t;
          for (int i = 0;i<3;i++) t[i]=r.gravity()(i);
          return t;
        },
        [](RBSystem& r, std::array<double,3> t){
          for (int i = 0;i<3;i++) r.gravity()(i)=t[i];
        })
      .def("beams", &RBSystem::beams)
      .def("springs", &RBSystem::springs)
      .def("saveState", &RBSystem::saveState)
      .def("reset", &RBSystem::reset)
      .def("connectorPos", [](RBSystem &r, Connector c){ auto v= r.connectorPos(c); return py::make_tuple(v(0),v(1),v(2));});
  
}
