// Copyright (c) 2026 Mikhail S. Kalenkov <mikhail.kalenkov@gmail.com>
// Licensed under the MIT License
#define BOOST_TEST_MODULE josephson test
#include "AlphaCoefficients.hpp"
#include "AndreevAmplitude.hpp"
#include "Config.hpp"
#include "EnergyGrid.hpp"
#include "MatrixConstants.hpp"
#include "ScatteringMatrix.hpp"
#include "Yfunction.hpp"
#include "misc.hpp"
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/program_options.hpp>
#include <complex>
#include <random>
#include <boost/dll/runtime_symbol_info.hpp>
#include <filesystem>
#include <fstream>

using namespace std::complex_literals;
using namespace std;
namespace po = boost::program_options;

double multiplier = 1.0e10;
bool fast_mode = false;
std::uint32_t seed = 0;
std::mt19937 rng;
size_t random_test_count = 100;

struct ArgsParser
{
    ArgsParser()
    {
        auto &mst = boost::unit_test::framework::master_test_suite();
        std::vector<std::string> args;
        for (int i = 1; i < mst.argc; ++i)
        {
            args.emplace_back(mst.argv[i]);
        }
        po::options_description desc("options");
        desc.add_options()                                     //
            ("fast", po::bool_switch(&fast_mode), "fast mode") //
            ("seed", po::value<std::uint32_t>(&seed), "random seed");

        po::variables_map vm;
        po::store(po::command_line_parser(args).options(desc).run(), vm);
        po::notify(vm);
        if (!vm.count("seed"))
        {
            seed = std::random_device{}();
        }
        rng = std::mt19937(seed);
        random_test_count = fast_mode ? 10 : 100;
    }
};

BOOST_GLOBAL_FIXTURE(ArgsParser);

BOOST_AUTO_TEST_CASE(test_matrix_commutation_relations)
{
    for (auto &m : {MatrixConstants::O3T3, MatrixConstants::O1, MatrixConstants::S1, MatrixConstants::S3, MatrixConstants::T1, MatrixConstants::I})
    {
        BOOST_TEST((m * m - MatrixConstants::I).norm() == 0);
    }
    BOOST_TEST((MatrixConstants::O3T3 * MatrixConstants::O1 + MatrixConstants::O1 * MatrixConstants::O3T3).norm() == 0);
    BOOST_TEST((MatrixConstants::O3T3 * MatrixConstants::S3 - MatrixConstants::S3 * MatrixConstants::O3T3).norm() == 0);
    BOOST_TEST((MatrixConstants::O3T3 * MatrixConstants::T1 + MatrixConstants::T1 * MatrixConstants::O3T3).norm() == 0);
    BOOST_TEST((MatrixConstants::O3T3 * MatrixConstants::S1 - MatrixConstants::S1 * MatrixConstants::O3T3).norm() == 0);
    BOOST_TEST((MatrixConstants::O1 * MatrixConstants::S3 - MatrixConstants::S3 * MatrixConstants::O1).norm() == 0);
    BOOST_TEST((MatrixConstants::O1 * MatrixConstants::T1 - MatrixConstants::T1 * MatrixConstants::O1).norm() == 0);
    BOOST_TEST((MatrixConstants::O1 * MatrixConstants::S1 - MatrixConstants::S1 * MatrixConstants::O1).norm() == 0);
    BOOST_TEST((MatrixConstants::T1 * MatrixConstants::S3 - MatrixConstants::S3 * MatrixConstants::T1).norm() == 0);
    BOOST_TEST((MatrixConstants::S1 * MatrixConstants::S3 + MatrixConstants::S3 * MatrixConstants::S1).norm() == 0);

    BOOST_TEST((MatrixConstants::A * MatrixConstants::S3 - MatrixConstants::S3 * MatrixConstants::A).norm() == 0);
    BOOST_TEST((MatrixConstants::A * MatrixConstants::A.adjoint() - MatrixConstants::I).norm() == 0);
    BOOST_TEST((MatrixConstants::Pp + MatrixConstants::Pm - MatrixConstants::I).norm() == 0);
    BOOST_TEST((MatrixConstants::Pp - MatrixConstants::Pm - MatrixConstants::O3T3).norm() == 0);
}

BOOST_AUTO_TEST_CASE(test_energy_grid)
{
    std::uniform_real_distribution<double> dist_energy_cutoff(2.0, 20.0);
    std::uniform_real_distribution<double> dist_log_voltage(std::log(1.0e-2), std::log(2.0));
    std::uniform_real_distribution<double> dist_log_step_size(std::log(1.0e-2), std::log(0.1));
    for (size_t c = 0; c < random_test_count; ++c)
    {
        double energy_cutoff = dist_energy_cutoff(rng);
        double voltage = std::exp(dist_log_voltage(rng));
        double energy_step_size = std::exp(dist_log_step_size(rng));
        EnergyGrid E(energy_cutoff, voltage, energy_step_size);
        auto step = E[1] - E[0];
        BOOST_TEST(std::abs(E[1] - E[0] - E.step_size()) * multiplier <= E.step_size());
        BOOST_TEST(std::abs(E.step_size() * E.steps_per_voltage() - E.voltage()) * multiplier <= 1.0);
        BOOST_TEST(std::abs(E.step_size() * E.steps_per_half_voltage() - E.voltage() / 2.0) * multiplier <= 1.0);
        BOOST_TEST(std::abs(E.energy_cutoff() - E.grid_max_energy()) >= voltage);
        for (size_t i = 0; i < E.size() - 1; i++)
        {
            BOOST_TEST(std::abs(E[i + 1] - E[i] - E.step_size()) * multiplier <= 1.0);
            BOOST_TEST(std::abs(E[i] + E[E.size() - 1 - i]) * multiplier <= 1.0);
            BOOST_TEST(std::abs(E[i] - E.at(i)) * multiplier <= 1.0);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_scattering_matrix)
{
    std::uniform_real_distribution<double> dist_reflection_probability_up(0.0, 1.0);
    std::uniform_real_distribution<double> dist_reflection_probability_down(0.0, 1.0);
    std::uniform_real_distribution<double> dist_spin_mixing_angle(0.0, 2.0 * M_PI);
    for (size_t c = 0; c < random_test_count; ++c)
    {
        double reflection_probability_up = dist_reflection_probability_up(rng);
        double reflection_probability_down = dist_reflection_probability_down(rng);
        double spin_mixing_angle = dist_spin_mixing_angle(rng);
        ScatteringMatrix S1(reflection_probability_up, reflection_probability_down, spin_mixing_angle);
        ScatteringMatrix S2(reflection_probability_down, reflection_probability_up, -spin_mixing_angle);
        BOOST_TEST((S1.Sn() * S1.Sn().adjoint() - Eigen::Matrix<std::complex<double>, 8, 8>::Identity()).norm() * multiplier <= 1.0);
        BOOST_TEST((S2.Sn() * S2.Sn().adjoint() - Eigen::Matrix<std::complex<double>, 8, 8>::Identity()).norm() * multiplier <= 1.0);
        BOOST_TEST((S1.Sn() - S1.Sn().transpose()).norm() * multiplier <= 1.0);
        BOOST_TEST((S2.Sn() - S2.Sn().transpose()).norm() * multiplier <= 1.0);
        BOOST_TEST((S1.Sn() * MatrixConstants::S3 - MatrixConstants::S3 * S1.Sn()).norm() * multiplier <= 1.0);
        BOOST_TEST((S2.Sn() * MatrixConstants::S3 - MatrixConstants::S3 * S2.Sn()).norm() * multiplier <= 1.0);
        BOOST_TEST((S1.Sn() - MatrixConstants::S1 * S2.Sn() * MatrixConstants::S1).norm() * multiplier <= 1.0);

        BOOST_TEST(std::abs(S1.reflection_probability_up() - S1.reflection_amplitude_up() * S1.reflection_amplitude_up()) * multiplier <= 1.0);
        BOOST_TEST(std::abs(S2.reflection_probability_up() - S2.reflection_amplitude_up() * S2.reflection_amplitude_up()) * multiplier <= 1.0);
        BOOST_TEST(std::abs(S1.reflection_probability_down() - S1.reflection_amplitude_down() * S1.reflection_amplitude_down()) * multiplier <= 1.0);
        BOOST_TEST(std::abs(S2.reflection_probability_down() - S2.reflection_amplitude_down() * S2.reflection_amplitude_down()) * multiplier <= 1.0);
        BOOST_TEST(std::abs(S1.transmission_probability_up() - S1.transmission_amplitude_up() * S1.transmission_amplitude_up()) * multiplier <= 1.0);
        BOOST_TEST(std::abs(S2.transmission_probability_up() - S2.transmission_amplitude_up() * S2.transmission_amplitude_up()) * multiplier <= 1.0);
        BOOST_TEST(std::abs(S1.transmission_probability_down() - S1.transmission_amplitude_down() * S1.transmission_amplitude_down()) * multiplier <= 1.0);
        BOOST_TEST(std::abs(S2.transmission_probability_down() - S2.transmission_amplitude_down() * S2.transmission_amplitude_down()) * multiplier <= 1.0);

        BOOST_TEST(std::abs(S1.reflection_probability_up() + S1.transmission_probability_up() - 1.0) * multiplier <= 1.0);
        BOOST_TEST(std::abs(S2.reflection_probability_up() + S2.transmission_probability_up() - 1.0) * multiplier <= 1.0);
        BOOST_TEST(std::abs(S1.reflection_probability_down() + S1.transmission_probability_down() - 1.0) * multiplier <= 1.0);
        BOOST_TEST(std::abs(S2.reflection_probability_down() + S2.transmission_probability_down() - 1.0) * multiplier <= 1.0);
    }
}

BOOST_AUTO_TEST_CASE(test_andreev_amplitude)
{
    std::uniform_real_distribution<double> dist_energy_cutoff(2.0, 20.0);
    std::uniform_real_distribution<double> dist_log_voltage(std::log(1.0e-2), std::log(2.0));
    std::uniform_real_distribution<double> dist_log_step_size(std::log(1.0e-2), std::log(0.1));
    std::uniform_real_distribution<double> dist_log_energy_decay(std::log(1.0e-8), std::log(2.0));
    std::uniform_real_distribution<double> dist_order_parameter(0.0, 2.0);
    std::uniform_real_distribution<double> dist_order_parameter_cutoff_smoothing(0.1, 2.0);
    for (size_t c = 0; c < random_test_count; ++c)
    {
        double energy_cutoff = dist_energy_cutoff(rng);
        double voltage = std::exp(dist_log_voltage(rng));
        double energy_step_size = std::exp(dist_log_step_size(rng));
        double energy_decay = std::exp(dist_log_energy_decay(rng));
        double order_parameter = dist_order_parameter(rng);
        double order_parameter_cutoff_smoothing = dist_order_parameter_cutoff_smoothing(rng);
        EnergyGrid E(energy_cutoff, voltage, energy_step_size);
        AndreevAmplitude a(E, energy_decay, order_parameter, order_parameter_cutoff_smoothing);

        for (int i = 0; i < E.size() / 2; i++)
        {
            BOOST_TEST(a.order_parameter(i) <= a.order_parameter(i + 1));
            BOOST_TEST(a.order_parameter(E.size() - 1 - i) <= a.order_parameter(E.size() - 1 - i - 1));
        };

        double derivative1 = 0;
        double derivative2 = 0;
        double max_derivative = 0;
        for (int i = 0; i < E.size() - 1; i++)
        {
            derivative1 = derivative2;
            if (E[i + 1] < -E.energy_cutoff())
            {
                derivative2 = 0;
            }
            if (E[i + 1] >= -E.energy_cutoff() and E[i + 1] < -E.energy_cutoff() + a.order_parameter_cutoff_smoothing())
            {
                derivative2 = (a.order_parameter() / a.order_parameter_cutoff_smoothing()) * sigmoid_derivative((E[i + 1] + E.energy_cutoff()) / a.order_parameter_cutoff_smoothing());
            }
            if (E[i + 1] >= -E.energy_cutoff() + a.order_parameter_cutoff_smoothing() and E[i + 1] < E.energy_cutoff() - a.order_parameter_cutoff_smoothing())
            {
                derivative2 = 0;
            }
            if (E[i + 1] >= E.energy_cutoff() - a.order_parameter_cutoff_smoothing() and E[i + 1] < E.energy_cutoff())
            {
                derivative2 = -(a.order_parameter() / a.order_parameter_cutoff_smoothing()) * sigmoid_derivative((E.energy_cutoff() - E[i + 1]) / a.order_parameter_cutoff_smoothing());
            }
            if (E[i + 1] >= E.energy_cutoff())
            {
                derivative2 = 0;
            }
            if (
                (E[i] < -E.energy_cutoff() + a.order_parameter_cutoff_smoothing() / 2.0 and E[i + 1] > -E.energy_cutoff() + a.order_parameter_cutoff_smoothing() / 2.0) or (E[i] < E.energy_cutoff() - a.order_parameter_cutoff_smoothing() / 2.0 and E[i + 1] > E.energy_cutoff() - a.order_parameter_cutoff_smoothing() / 2.0))
            {
                max_derivative = 1.5 * (a.order_parameter() / a.order_parameter_cutoff_smoothing());
            }
            else
            {
                max_derivative = std::max(std::abs(derivative1), std::abs(derivative2));
            }
            BOOST_TEST(std::abs(a.order_parameter(i) - a.order_parameter(i + 1)) <= (E[i + 1] - E[i]) * max_derivative);
        }
        for (int i = 0; i < E.size(); i++)
        {
            BOOST_TEST(std::abs(a.order_parameter(i) - a.order_parameter(E.size() - 1 - i)) * multiplier <= 1.0);
            BOOST_TEST(std::abs(a[i]) <= 1.0);
            BOOST_TEST(std::abs(a.order_parameter(i) * a[i] * a[i] + 2.0 * (E[i] + 1i * a.energy_decay()) * a[i] + a.order_parameter(i)) * multiplier <= 1.0);
            if (a.order_parameter(i) == 0)
            {
                BOOST_TEST(a[i].imag() >= 0);
            }
            else
            {
                BOOST_TEST(a[i].imag() > 0);
            }
        }
        int N = 100 * E.steps_per_voltage();
        for (int i = 0 - N; i < E.size() + N; i++)
        {
            if (i >= 0 and i < E.size())
            {
                BOOST_TEST(std::abs(a.order_parameter(i) - a.order_parameter_at(i)) * multiplier <= 1.0);
                BOOST_TEST(std::abs(a[i] - a.at(i)) * multiplier <= 1.0);
            }
            else
            {
                BOOST_TEST(std::abs(a.order_parameter_at(i)) * multiplier <= 1.0);
                BOOST_TEST(std::abs(a.at(i)) * multiplier <= 1.0);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_alpha_coeffcients)
{
    std::uniform_real_distribution<double> dist_energy_cutoff(2.0, 20.0);
    std::uniform_real_distribution<double> dist_log_voltage(std::log(1.0e-2), std::log(2.0));
    std::uniform_real_distribution<double> dist_log_step_size(std::log(1.0e-2), std::log(0.1));
    std::uniform_real_distribution<double> dist_log_energy_decay(std::log(1.0e-8), std::log(2.0));
    std::uniform_real_distribution<double> dist_order_parameter(0.0, 2.0);
    std::uniform_real_distribution<double> dist_order_parameter_cutoff_smoothing(0.1, 2.0);

    std::uniform_real_distribution<double> dist_reflection_probability_up(0.0, 1.0);
    std::uniform_real_distribution<double> dist_reflection_probability_down(0.0, 1.0);
    std::uniform_real_distribution<double> dist_spin_mixing_angle(0.0, 2.0 * M_PI);
    for (size_t c = 0; c < random_test_count; ++c)
    {
        double energy_cutoff = dist_energy_cutoff(rng);
        double voltage = std::exp(dist_log_voltage(rng));
        double energy_step_size = std::exp(dist_log_step_size(rng));
        double energy_decay = std::exp(dist_log_energy_decay(rng));
        double order_parameter = dist_order_parameter(rng);
        double order_parameter_cutoff_smoothing = dist_order_parameter_cutoff_smoothing(rng);
        double reflection_probability_up = dist_reflection_probability_up(rng);
        double reflection_probability_down = dist_reflection_probability_down(rng);
        double spin_mixing_angle = dist_spin_mixing_angle(rng);
        EnergyGrid E(energy_cutoff, voltage, energy_step_size);
        AndreevAmplitude a(E, energy_decay, order_parameter, order_parameter_cutoff_smoothing);
        ScatteringMatrix S(reflection_probability_up, reflection_probability_down, spin_mixing_angle);
        AlphaCoefficients alpha(E, S, a);
        Eigen::Matrix<std::complex<double>, 8, 8> diff_alpha_p;
        Eigen::Matrix<std::complex<double>, 8, 8> diff_alpha_m;
        Eigen::Matrix<std::complex<double>, 8, 8> diff_alpha_0;
        int N = 100 * E.steps_per_voltage();
        for (int i = 0; i < E.size(); i++)
        {
            if (i + E.steps_per_voltage() < E.size())
            {
                diff_alpha_p = alpha.p(i) * alpha.p(i + E.steps_per_voltage()) * a[i + E.steps_per_voltage()] * MatrixConstants::Pm + alpha.p(i) * S.W() + a[i] * MatrixConstants::Pp;
            }
            else
            {
                diff_alpha_p = alpha.p(i) * S.W() + a[i] * MatrixConstants::Pp;
            }
            BOOST_TEST(diff_alpha_p.norm() * multiplier < 1.0);
        }
        for (int i = 0 - N; i < E.size() + N; i++)
        {
            diff_alpha_p = alpha.p_at(i) * alpha.p_at(i + E.steps_per_voltage()) * a.at(i + E.steps_per_voltage()) * MatrixConstants::Pm + alpha.p_at(i) * S.W() + a.at(i) * MatrixConstants::Pp;

            BOOST_TEST(diff_alpha_p.norm() * multiplier < 1.0);
        }
        for (int i = 0; i < E.size(); i++)
        {
            if (i - E.steps_per_voltage() >= 0)
            {
                diff_alpha_m = alpha.m(i) * alpha.m(i - E.steps_per_voltage()) * a[i - E.steps_per_voltage()] * MatrixConstants::Pp + alpha.m(i) * S.W() + a[i] * MatrixConstants::Pm;
            }
            else
            {
                diff_alpha_m = alpha.m(i) * S.W() + a[i] * MatrixConstants::Pm;
            }
            BOOST_TEST(diff_alpha_m.norm() * multiplier < 1.0);
        }
        for (int i = 0 - N; i < E.size() + N; i++)
        {
            diff_alpha_m = alpha.m_at(i) * alpha.m_at(i - E.steps_per_voltage()) * a.at(i - E.steps_per_voltage()) * MatrixConstants::Pp + alpha.m_at(i) * S.W() + a.at(i) * MatrixConstants::Pm;

            BOOST_TEST(diff_alpha_m.norm() * multiplier < 1.0);
        }
        for (int i = 0; i < E.size(); i++)
        {
            if (i - E.steps_per_half_voltage() >= 0 and i + E.steps_per_half_voltage() < E.size())
            {
                diff_alpha_0 = alpha.z(i) * (S.W()                                                                                               //
                                             + alpha.m(i - E.steps_per_half_voltage()) * a[i - E.steps_per_half_voltage()] * MatrixConstants::Pp //
                                             + alpha.p(i + E.steps_per_half_voltage()) * a[i + E.steps_per_half_voltage()] * MatrixConstants::Pm) -
                               MatrixConstants::I;
            }
            if (i - E.steps_per_half_voltage() < 0)
            {
                diff_alpha_0 = alpha.z(i) * (S.W() //
                                             + alpha.p(i + E.steps_per_half_voltage()) * a[i + E.steps_per_half_voltage()] * MatrixConstants::Pm) -
                               MatrixConstants::I;
            }
            if (i + E.steps_per_half_voltage() >= E.size())
            {
                diff_alpha_0 = alpha.z(i) * (S.W() //
                                             + alpha.m(i - E.steps_per_half_voltage()) * a[i - E.steps_per_half_voltage()] * MatrixConstants::Pp) -
                               MatrixConstants::I;
            }
            BOOST_TEST(diff_alpha_0.norm() * multiplier < 1.0);
        }
        for (int i = 0 - N; i < E.size() + N; i++)
        {
            diff_alpha_0 = alpha.z_at(i) * (S.W()                                                                                                     //
                                            + alpha.m_at(i - E.steps_per_half_voltage()) * a.at(i - E.steps_per_half_voltage()) * MatrixConstants::Pp //
                                            + alpha.p_at(i + E.steps_per_half_voltage()) * a.at(i + E.steps_per_half_voltage()) * MatrixConstants::Pm) -
                           MatrixConstants::I;
            BOOST_TEST(diff_alpha_0.norm() * multiplier < 1.0);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_y_function)
{
    std::uniform_real_distribution<double> dist_energy_cutoff(2.0, 20.0);
    std::uniform_real_distribution<double> dist_log_voltage(std::log(1.0e-2), std::log(2.0));
    std::uniform_real_distribution<double> dist_log_step_size(std::log(1.0e-2), std::log(0.1));
    std::uniform_real_distribution<double> dist_log_energy_decay(std::log(1.0e-8), std::log(2.0));
    std::uniform_real_distribution<double> dist_order_parameter(0.0, 2.0);
    std::uniform_real_distribution<double> dist_order_parameter_cutoff_smoothing(0.1, 2.0);

    std::uniform_real_distribution<double> dist_reflection_probability_up(0.0, 1.0);
    std::uniform_real_distribution<double> dist_reflection_probability_down(0.0, 1.0);
    std::uniform_real_distribution<double> dist_spin_mixing_angle(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<double> dist_superconducting_phase(0.0, 2.0 * M_PI);
    for (size_t c = 0; c < random_test_count; ++c)
    {
        double energy_cutoff = dist_energy_cutoff(rng);
        double voltage = std::exp(dist_log_voltage(rng));
        double energy_step_size = std::exp(dist_log_step_size(rng));
        double energy_decay = std::exp(dist_log_energy_decay(rng));
        double order_parameter = dist_order_parameter(rng);
        double order_parameter_cutoff_smoothing = dist_order_parameter_cutoff_smoothing(rng);
        double reflection_probability_up = dist_reflection_probability_up(rng);
        double reflection_probability_down = dist_reflection_probability_down(rng);
        double spin_mixing_angle = dist_spin_mixing_angle(rng);
        double superconducting_phase = dist_superconducting_phase(rng);
        EnergyGrid E(energy_cutoff, voltage, energy_step_size);
        AndreevAmplitude a(E, energy_decay, order_parameter, order_parameter_cutoff_smoothing);
        ScatteringMatrix S(reflection_probability_up, reflection_probability_down, spin_mixing_angle);
        AlphaCoefficients alpha(E, S, a);
        Yfunction Y(alpha, E, S, a, superconducting_phase);
        std::atomic<std::size_t> progress;
        Y.evaluateY(progress);
        Eigen::Matrix<std::complex<double>, 8, 8> diff_Y;
        for (int i = 2 * E.steps_per_half_voltage(); i < E.size() - 2 * E.steps_per_half_voltage(); i++)
        {
            diff_Y = Y[i] * S.W() + std::exp(-1.0i * superconducting_phase / 2.0) * Y[i + 2 * E.steps_per_half_voltage()] * a[i + E.steps_per_half_voltage()] * MatrixConstants::Pm + std::exp(1.0i * superconducting_phase / 2.0) * Y[i - 2 * E.steps_per_half_voltage()] * a[i - E.steps_per_half_voltage()] * MatrixConstants::Pp - MatrixConstants::I;
            BOOST_TEST(diff_Y.norm() * multiplier < 1.0);
        }
        int N = 100 * E.steps_per_voltage();
        for (int i = 0 - N; i < E.size() + N; i++)
        {
            diff_Y = Y.at(i) * S.W() + std::exp(-1.0i * superconducting_phase / 2.0) * Y.at(i + 2 * E.steps_per_half_voltage()) * a.at(i + E.steps_per_half_voltage()) * MatrixConstants::Pm + std::exp(1.0i * superconducting_phase / 2.0) * Y.at(i - 2 * E.steps_per_half_voltage()) * a.at(i - E.steps_per_half_voltage()) * MatrixConstants::Pp - MatrixConstants::I;
            BOOST_TEST(diff_Y.norm() * multiplier < 1.0);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_josephson_current)
{
    auto exe_dir = boost::dll::program_location().parent_path();
    auto datafile = exe_dir / "../test/data.csv";
    std::ifstream file(datafile);
    if (!file)
        throw std::runtime_error("Cannot open file: " + datafile.string());
    std::string line;
    // process header
    std::getline(file, line);
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ','))
    {
        fields.push_back(item);
    }
    // process the rest of the file
    while (std::getline(file, line))
    {
        std::vector<std::string> options;
        std::stringstream ss(line);
        std::string item;
        while (std::getline(ss, item, ','))
        {
            options.push_back(item);
        }
        Config config;
        config.R_up(std::stod(options[0]));
        config.R_down(std::stod(options[1]));
        config.spin_mixing_angle(std::stod(options[2]));
        config.energy_cutoff(std::stod(options[3]));
        config.step_size(std::stod(options[6]));
        config.eV(std::stod(options[4]));
        config.T(std::stod(options[5]));
        config.energy_broadening(std::stod(options[7]));
        config.order_parameter_cutoff_smoothing(std::stod(options[9]));
        if (fast_mode and (options[11] != "fast"))
        {
            continue;
        }
        EnergyGrid E(config.energy_cutoff(), config.eV(), config.step_size());
        ScatteringMatrix S(config.R_up(), config.R_down(), config.spin_mixing_angle());
        AndreevAmplitude a(E, config.energy_broadening(), 1.0, config.order_parameter_cutoff_smoothing());
        AlphaCoefficients alpha(E, S, a);
        std::vector J = evaluateJphi(alpha, E, S, a, config, {std::stod(options[8])});
        double Jsaved = std::stod(options[10]);
        BOOST_TEST(std::abs(J[0] - Jsaved) * multiplier < 1.0);
    }
}