#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

#include "rbs_fem_SR.h"
#include "rigid_body_helper_SR.h"
#include "eq_system_fem_SR.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(Transformation<>);
// PYBIND11_MAKE_OPAQUE(MassMatrix);

PYBIND11_MODULE(rigid_body_FEM_SR, rbdsr) {

    // adds the BLA bindings as submodule, accessible as rigid_body.bla.Matrix etc.
    // https://github.com/pybind/pybind11/discussions/4027
    auto m = rbdsr.def_submodule("bla", "basic linear algebra");
    #include "bind_bla_obj.h"



    // the main bindings:
    rbdsr.doc() = "rigid body FEM Shake and Rattle simulator";
       
    py::class_<Transformation<>>(rbdsr, "Transformation", py::module_local())
      .def(py::init<>())
      .def(py::init<Vector<double>>())
      .def("__str__", [](Transformation<> & t) {
        std::stringstream sstr;
        sstr << t;
        return sstr.str();
      })
      .def("setTranslation", py::overload_cast<double, double, double>(&Transformation<>::setTranslation))
      .def("setTranslation", py::overload_cast<VectorView<double>>(&Transformation<>::setTranslation))

      .def("setRotation", py::overload_cast<size_t, size_t, double>(&Transformation<>::setRotation))
      .def("setRotation", py::overload_cast<MatrixView<double>>(&Transformation<>::setRotation))

      .def("asTuple",[](Transformation<>& t){
        // *column-major* transformation matrix as in https://threejs.org/docs/#api/en/math/Matrix4
        // converts Schöberl-style ordering of Q to column-major ordering of a three.js transformation matrix:
        return py::make_tuple(t.q_(3), t.q_(6), t.q_(9), 0, t.q_(4), t.q_(7), t.q_(10), 0, t.q_(5), t.q_(8), t.q_(11), 0, t.q_(0), t.q_(1), t.q_(2), 1);
      });


    py::class_<Connector>(rbdsr,"Connector", py::module_local())
      .def_property("pos",[](Connector& c){return py::make_tuple(c.Pos(0),c.Pos(1),c.Pos(2));},
                          [](Connector& c, std::array<double,3> t){c.Pos(0) = t[0];c.Pos(1)=t[1];c.Pos(2)=t[2];})
      .def_property_readonly("body_index",[](Connector& c){return c.Body_index();})
      .def_property_readonly("type",[](Connector& c){return c.Type() == ConnectorType::mass ? 0 : 1 ;});

    py::class_<Beam>(rbdsr,"Beam", py::module_local())
      .def(py::init<>([](Connector a, Connector b){return Beam{a,b};}))
      .def_property_readonly("length", [](Beam& b){return b.Length();})
      .def_property_readonly("connectorA", [](Beam& b){return b.Connector_a();})
      .def_property_readonly("connectorB",[](Beam& b){return b.Connector_b();});
    
    py::class_<Spring>(rbdsr,"Spring", py::module_local())
      .def(py::init<>([](Connector a, Connector b, double length, double stiffness){return Spring{length,stiffness,a,b};})) // stiffness should be positive! (-k)
      .def_property_readonly("length", [](Spring& b){return b.Length();})
      .def_property_readonly("stiffness", [](Spring& b){return b.Stiffness();}) // stiffness should be positive! (-k)
      .def_property_readonly("connectorA", [](Spring& b){return b.Connector_a();})
      .def_property_readonly("connectorB",[](Spring& b){return b.Connector_b();});

    py::class_<RigidBody_FEM> (rbdsr, "RigidBody_FEM_SR")
      .def(py::init<>())
      .def_property("transformation", &RigidBody_FEM::getQ, &RigidBody_FEM::setQ)
      // .def_property("phat", &RigidBody_FEM::getPhat,&RigidBody_FEM::setPhat)
      .def_property("momentumTrans",
        [](RigidBody_FEM& r){
          return py::make_tuple(r.phat()(0), r.phat()(1), r.phat()(2));
        },
        [](RigidBody_FEM& r, std::array<double,3> vals){
          for (int i=0; i < 3; i++)
            r.phat()(i) = vals[i];
        })
      .def_property("momentumRot",
        [](RigidBody_FEM& r){
          return py::make_tuple(r.phat()(3), r.phat()(4), r.phat()(5));
        },
        [](RigidBody_FEM& r, std::array<double,3> vals){
          for (int i=0; i < 3; i++)
            r.phat()(i + 3) = vals[i];
        })
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
      .def("reset", &RigidBody_FEM::reset);
      // .def("setPhat", &RigidBody_FEM::setPhat_v)
      //.def(py::pickle(
      //    [](RigidBody_FEM& rbdsr){ // __getstate__
      //        return rbdsr.to_pickle();
      //    },
      //    [](py::tuple data){ // __setstate__
      //      RigidBody_FEM rbdsr;
      //      rbdsr.load_pickle(data);
      //      return rbdsr;
      //    }));
      
    //rbdsr.def("mass_matrix_from_inertia", &mass_matrix_from_inertia, "generates the a mass matrix from given inertia, center and mass",
            //py::arg("inertia_matrix"), py::arg("center_of_mass"), py::arg("mass"));


    py::class_<RBS_FEM> (rbdsr, "RBS_FEM_SR")
      .def(py::init<>())
      .def("saveState", &RBS_FEM::saveState)
      .def("add", py::overload_cast<RigidBody_FEM&>(&RBS_FEM::add))
      .def("add", py::overload_cast<Beam&>(&RBS_FEM::add))
      .def("add", py::overload_cast<Spring&>(&RBS_FEM::add))
      .def("addFix",&RBS_FEM::addFix)
      .def("bodies", py::overload_cast<>(&RBS_FEM::Bodies), py::return_value_policy::reference_internal)
      .def("bodies", py::overload_cast<size_t>(&RBS_FEM::Bodies), py::return_value_policy::reference_internal)
      .def_property("gravity",
        [](RBS_FEM& r){
          std::array<double,3> t;
          for (int i = 0;i<3;i++) t[i]=r.Gravity()(i);
          return t;
        },
        [](RBS_FEM& r, std::array<double,3> t){
          for (int i = 0;i<3;i++) r.Gravity()(i)=t[i];
        })
      .def("beams", py::overload_cast<>(&RBS_FEM::Beams))
      .def("beams", py::overload_cast<size_t>(&RBS_FEM::Beams))
      .def("springs", py::overload_cast<>(&RBS_FEM::Springs))
      .def("springs", py::overload_cast<size_t>(&RBS_FEM::Springs))
      .def("saveState", &RBS_FEM::saveState)
      .def("reset", &RBS_FEM::reset)
      .def("connectorPos", [](RBS_FEM &r, Connector c){auto v = r.connectorPos(c); return py::make_tuple(v(0),v(1),v(2));});


    rbdsr.def("simulate",[](RBS_FEM& rbs, double tend, double steps) {simulate(rbs, tend, steps);});
}
