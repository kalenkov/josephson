// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Config.hpp"
#include "JosephsonCurrent.hpp"

namespace py = pybind11;

Config to_cpp_config(const py::dict &pycfg)
{
    Config cfg;
    cfg.R_up(pycfg["R_up"].cast<double>());
    cfg.R_down(pycfg["R_down"].cast<double>());
    cfg.spin_mixing_angle(pycfg["spin_mixing_angle"].cast<double>());
    cfg.energy_cutoff(pycfg["energy_cutoff"].cast<double>());
    cfg.step_size(pycfg["step_size"].cast<double>());
    cfg.eV(pycfg["eV"].cast<double>());
    cfg.T(pycfg["T"].cast<double>());
    cfg.energy_broadening(pycfg["energy_broadening"].cast<double>());
    cfg.order_parameter_cutoff_smoothing(pycfg["order_parameter_cutoff_smoothing"].cast<double>());
    cfg.show_progress(pycfg["show_progress"].cast<bool>());
    return cfg;
}

PYBIND11_MODULE(josephson, m)
{
    m.doc() = "Josephson current module";

    py::class_<RuntimeInfo>(m, "RuntimeInfo")
        .def_readonly("step_size", &RuntimeInfo::step_size)
        .def_readonly("grid_max_energy", &RuntimeInfo::grid_max_energy)
        .def_readonly("grid_size", &RuntimeInfo::grid_size)
        .def("to_dict", [](const RuntimeInfo &r)
        {
            py::dict d;
            d["step_size"] = r.step_size;
            d["grid_max_energy"] = r.grid_max_energy;
            d["grid_size"] = r.grid_size;
            return d;
        }); 

    py::class_<JosephsonCurrent>(m, "JosephsonCurrent")
        .def(py::init([](const py::dict &pycfg)
        {
            return JosephsonCurrent(to_cpp_config(pycfg));
        }))
        .def("evaluate", py::overload_cast<double>(&JosephsonCurrent::evaluate, py::const_))
        .def("evaluate", py::overload_cast<const std::vector<double> &>(&JosephsonCurrent::evaluate, py::const_))
        .def("runtime_info", &JosephsonCurrent::runtime_info);
}