// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#include "AlphaCoefficients.hpp"
#include "AndreevAmplitude.hpp"
#include "Config.hpp"
#include "EnergyGrid.hpp"
#include "ScatteringMatrix.hpp"
#include "Yfunction.hpp"
#include "RuntimeInfo.hpp"
#include <format>
#include <omp.h>
#include <ostream>
#include <thread>

void print_config(std::ostream& out, const Config& config, const RuntimeInfo& rt_info, const double& phi_min, const double& phi_max, const int& n_phi, const std::string& output_filename, const std::string& cmd_line)
{
    out << std::format("# {:<12} = {}\n", "Command line", cmd_line)
        << "# Kernel configuration:\n"
        << std::format("# {:<12} = {}\n", "R-up", config.R_up())
        << std::format("# {:<12} = {}\n", "R-down", config.R_down())
        << std::format("# {:<12} = {}\n", "spin-mixing-angle", config.spin_mixing_angle())
        << std::format("# {:<12} = {}\n", "energy-cutoff", config.energy_cutoff())
        << std::format("# {:<12} = {}\n", "step-size", config.step_size())
        << std::format("# {:<12} = {}\n", "eV", config.eV())
        << std::format("# {:<12} = {}\n", "T", config.T())
        << std::format("# {:<12} = {}\n", "energy-broadening", config.energy_broadening())
        << std::format("# {:<12} = {}\n", "order-parameter-cutoff-smoothing", config.order_parameter_cutoff_smoothing())
        << std::format("# {:<12} = {}\n", "show-progress", config.show_progress())
        << "# Calculated parameters:\n"
        << std::format("# {:<12} = {}\n", "maximal grid energy", rt_info.grid_max_energy)
        << std::format("# {:<12} = {}\n", "calculated step size", rt_info.step_size)
        << std::format("# {:<12} = {}\n", "grid size", rt_info.grid_size)
        << "# Program configuration:\n"
        << std::format("# {:<12} = {}\n", "phi-min", phi_min)
        << std::format("# {:<12} = {}\n", "phi-max", phi_max)
        << std::format("# {:<12} = {}\n", "n-phi", n_phi)
        << std::format("# {:<12} = {}\n", "filename", output_filename)
        << std::endl;
}

std::vector<double> evaluateJphi(const AlphaCoefficients& alpha,
    const EnergyGrid& E, const ScatteringMatrix& S,
    const AndreevAmplitude& a, const Config& config, std::vector<double> phi)
{
    const int N = phi.size();
    const int steps_per_eval = (E.size() + (E.num_voltage_steps() - 1) * (2 * E.size() - E.num_voltage_steps() * E.steps_per_voltage())) / 1000;
    const std::size_t total_steps = static_cast<std::size_t>(N) * steps_per_eval;

    std::atomic<std::size_t> completed_steps { 0 };
    std::atomic<bool> done { false };

    std::vector<double> J(N);

    std::thread progress_thread([&]() {
        auto start_time = std::chrono::steady_clock::now();
        while (!done) {
            auto now = std::chrono::steady_clock::now();
            std::size_t done_steps = completed_steps.load();
            double percent = 100.0 * done_steps / total_steps;

            double elapsed_seconds = std::chrono::duration<double>(now - start_time).count();
            double estimated_total_seconds = (done_steps > 0)
                ? (elapsed_seconds * total_steps / done_steps)
                : 0.0;
            double remaining_seconds = estimated_total_seconds - elapsed_seconds;

            int elapsed_min = static_cast<int>(elapsed_seconds) / 60;
            int elapsed_sec = static_cast<int>(elapsed_seconds) % 60;
            int remaining_min = static_cast<int>(remaining_seconds) / 60;
            int remaining_sec = static_cast<int>(remaining_seconds) % 60;
            if (config.show_progress())
            {
                std::cout << std::format("Progress: {:5.1f}% ({}/{}) | Elapsed: {:02d}:{:02d} | ETA: {:02d}:{:02d}\r",
                                         percent, done_steps, total_steps,
                                         elapsed_min, elapsed_sec,
                                         remaining_min, remaining_sec);
                std::cout.flush();
            }

            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    });

#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        Yfunction Y(alpha, E, S, a, phi[i]);
        Y.evaluateY(completed_steps);
        Y.evaluateQ(config.T());
        J[i] = Y.ElectricCurrent();
    }

    done = true;
    progress_thread.join();

    if (config.show_progress())
    {
        std::cout << std::format("Progress: 100.0% ({}/{})\n", completed_steps.load(), total_steps);
        std::cout.flush();
        std::cout << std::endl;
    }
    return J;
}

void printJphi(std::ostream& out, const std::vector<double>& phi, const std::vector<double>& J)
{
    out << "# phi J" << std::endl;
    for (int i = 0; i < phi.size(); i++) {
        out << std::format("{} {}", phi[i], J[i]) << std::endl;
    }
}
