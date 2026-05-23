// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#include "Config.hpp"
#include "misc.hpp"
#include "JosephsonCurrent.hpp"
#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>

using namespace std::complex_literals;
namespace po = boost::program_options;

int main(int argc, char** argv)
{
    Config cfg;
    std::ostringstream cmd;
    double R_up = cfg.R_up(), R_down = cfg.R_down(), spin_mixing_angle = cfg.spin_mixing_angle();
    double energy_cutoff = cfg.energy_cutoff(), eV = cfg.eV(), step_size = cfg.step_size();
    double energy_broadening = cfg.energy_broadening(), order_parameter_cutoff_smoothing = cfg.order_parameter_cutoff_smoothing();
    double T = cfg.T();
    double phi_min=0, phi_max=2.0 * M_PI;
    bool show_progress = false;
    int n_phi = 100;
    std::string output_filename("output.txt");
    for (int i = 0; i < argc; ++i)
    {
        if (i > 0)
            cmd << " ";
        cmd << argv[i];
    }

    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "Show help message") //
        ("R-up", po::value<double>(&R_up)->default_value(0.1), "Spin-up reflection probability") //
        ("R-down", po::value<double>(&R_down)->default_value(0.1), "Spin-down reflection probability") //
        ("spin-mixing-angle", po::value<double>(&spin_mixing_angle)->default_value(0.0), "Spin-mixing angle") //
        ("energy-cutoff", po::value<double>(&energy_cutoff)->default_value(20), "Energy cutoff") //
        ("step-size", po::value<double>(&step_size)->default_value(0.001), "Energy mesh step size") //
        ("eV", po::value<double>(&eV)->default_value(0.1), "Applied voltage") //
        ("T", po::value<double>(&T)->default_value(0.1), "Temperature") //
        ("energy-broadening", po::value<double>(&energy_broadening)->default_value(1e-07), "Decay parameter") //
        ("order-parameter-cutoff-smoothing", po::value<double>(&order_parameter_cutoff_smoothing)->default_value(2), "Smoothing width of the order-parameter cutoff at high energies") //
        ("phi-min", po::value<double>(&phi_min)->default_value(0.0), "Minimal phase") //
        ("phi-max", po::value<double>(&phi_max)->default_value(6.283185307179586), "Maximal phase") //
        ("n-phi", po::value<int>(&n_phi)->default_value(100), "Number of phase points") //
        ("output,o", po::value<std::string>(&output_filename)->default_value("output.txt"), "Output filename") //
        ("show-progress", po::bool_switch(&show_progress), "Show computation progress"); //

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        exit(0);
    }
    cfg.R_up(R_up);
    cfg.R_down(R_down);
    cfg.spin_mixing_angle(spin_mixing_angle);
    cfg.energy_cutoff(energy_cutoff);
    cfg.step_size(step_size);
    cfg.eV(eV);
    cfg.T(T);
    cfg.energy_broadening(energy_broadening);
    cfg.order_parameter_cutoff_smoothing(order_parameter_cutoff_smoothing);
    cfg.show_progress(show_progress);

    JosephsonCurrent jc(cfg);

    std::ofstream file(output_filename);
    if (!file) {
        std::cerr << "Error: Cannot open file " << output_filename << std::endl;
    }
    
    print_config(std::cout, cfg, jc.runtime_info(), phi_min, phi_max, n_phi, output_filename, cmd.str());
    print_config(file, cfg, jc.runtime_info(), phi_min, phi_max, n_phi, output_filename, cmd.str());

    std::vector<double> phi(n_phi);
    if (n_phi == 1)
    {
        phi[0] = phi_min;
    }
    else
    {
        for (int i = 0; i < n_phi; i++)
        {
            phi[i] = phi_min+ (phi_max - phi_min) * i / (n_phi - 1.0);
        }
    }
    std::vector<double> J = jc.evaluate(phi);

    printJphi(std::cout, phi, J);
    printJphi(file, phi, J);

    return 0;
}
