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
        // converts Schöberl-style ordering of Q to column-major ordering of a three.js transformation matrix:
        return py::make_tuple(t.q_(1), t.q_(5), t.q_(9), 0, t.q_(2), t.q_(6), t.q_(10), 0, t.q_(3), t.q_(7), t.q_(11), 0, t.q_(0), t.q_(4), t.q_(8), 1);
      });

    py::class_<Connector>(rbd,"Connector")
      .def_property("pos",[](Connector& c){return py::make_tuple(c.pos(0),c.pos(1),c.pos(2));},
                          [](Connector& c, std::array<double,3> t){c.pos(0)=t[0];c.pos(1)=t[1];c.pos(2)=t[2];})
      .def_property_readonly("body_index",[](Connector& c){return c.body_index;})
      .def_property_readonly("type",[](Connector& c){return c.t == ConnectorType::mass ? 0 : 1 ;});
    
    py::class_<Spring>(rbd,"Spring")
      .def(py::init<>([](Connector a, Connector b, double length, double stiffness){return Spring{length,stiffness,a,b};})) // stiffness should be positive! (-k)
      .def_property_readonly("length", [](Spring& b){return b.length;})
      .def_property_readonly("stiffness", [](Spring& b){return b.length;}) // stiffness should be positive! (-k)
      .def_property_readonly("connectorA", [](Spring& b){return b.a;})
      .def_property_readonly("connectorB",[](Spring& b){return b.b;});

    py::class_<RigidBody_FEM> (rbd, "RigidBody_FEM", py::dynamic_attr())
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
      .def_property("vertices",
        [](RigidBody_FEM& r){
          return r.vertices();
        },
        [](RigidBody_FEM& r, py::list v){
          r.vertices()=v;
      })
      .def_property("normals",
        [](RigidBody_FEM& r){
          return r.normals();
        },
        [](RigidBody_FEM& r, py::list n){
          r.normals()=n;
      })
      //.def("setMass", &RigidBody::setMass)
      .def("recalcMassMatrix", &RigidBody_FEM::recalcMassMatrix)
      .def("saveState", &RigidBody_FEM::saveState)
      .def("reset", &RigidBody_FEM::reset) 
      .def("setPhat", &RigidBody_FEM::setPhat_v);
      
    //rbd.def("mass_matrix_from_inertia", &mass_matrix_from_inertia, "generates the a mass matrix from given inertia, center and mass",
            //py::arg("inertia_matrix"), py::arg("center_of_mass"), py::arg("mass"));


    py::class_<RBS_FEM> (rbd, "RBS_FEM")
      .def(py::init<>())
      .def("saveState", &RBS_FEM::saveState)
      .def("addBody",&RBS_FEM::addBody)
      //.def("addBeam",&RBS_FEM::addBeam)
      .def("addSpring",&RBS_FEM::addSpring)
      .def("addFix",&RBS_FEM::addFix)
      .def("bodies", &RBS_FEM::bodies)
      .def_property("gravity",
        [](RBS_FEM& r){
          std::array<double,3> t;
          for (int i = 0;i<3;i++) t[i]=r.gravity()(i);
          return t;
        },
        [](RBS_FEM& r, std::array<double,3> t){
          for (int i = 0;i<3;i++) r.gravity()(i)=t[i];
        })
      //.def("beams", &RBS_FEM::beams)
      .def("springs", &RBS_FEM::springs)
      .def("saveState", &RBS_FEM::saveState)
      .def("reset", &RBS_FEM::reset)
      .def("connectorPos", [](RBS_FEM &r, Connector c){auto v = r.connectorPos(c); return py::make_tuple(v(0),v(1),v(2));});


    rbd.def("simulate",[](RBS_FEM& rbs, double tend, double steps) {simulate(rbs, tend, steps);});
  
}
