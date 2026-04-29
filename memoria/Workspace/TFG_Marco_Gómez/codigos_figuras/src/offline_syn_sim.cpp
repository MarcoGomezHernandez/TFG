#include <iostream>
#include <fstream>
#include <yaml-cpp/yaml.h>
#include "aux/scaling.hpp"
#include "aux/utils.hpp"
#include <ChemicalSynapsis.h>

typedef RungeKutta4 Integrator;
typedef HindmarshRoseNeuron<Integrator> NeuronType;
typedef ChemicalSynapsis<NeuronType, NeuronType, Integrator, double> ChemicalSynapseType;

int main(int argc, char *argv[])
{
    if (argc < 12)
    {
        std::cerr << "Usage: " << argv[0] << " <csv_path> <column_idx> <csv_step> <start_time> <stabilization_time> <use_time> <observation_time> <syn_component> <subsampling> <params_yaml_path> <output_csv_path>\n";
        std::cerr << "syn_component: 0=i_fast, 1=i_slow, 2=both (i, i_fast, i_slow)\n";
        return 1;
    }

    const std::string csv_path = argv[1];
    const size_t column_idx = static_cast<size_t>(std::stoul(argv[2]));
    const double csv_step = std::stod(argv[3]);
    const double start_time = std::stod(argv[4]);
    const double stabilization_time = std::stod(argv[5]);
    const double use_time = std::stod(argv[6]);
    const double observation_time = std::stod(argv[7]);
    const size_t syn_component = static_cast<size_t>(std::stoul(argv[8]));
    const size_t subsampling = static_cast<size_t>(std::stoul(argv[9]));
    const std::string params_yaml_path = argv[10];
    const std::string output_csv_path = argv[11];

    YAML::Node p = YAML::LoadFile(params_yaml_path);

    std::optional<ScaledSigResult> scaled_result = scale_sig(
        csv_path, column_idx, csv_step, start_time, stabilization_time + use_time, observation_time,
        NumericIntegrator::RK4, NeuronModel::HINDMARSH_ROSE, false);

    if (!scaled_result)
    {
        std::cerr << "Error al escalar la señal.\n";
        return 1;
    }

    NeuronType model_neur = create_hindmarsh_rose<Integrator>(false);
    reset_state_hindmarsh_rose(model_neur);

    const bool use_i_fast = (syn_component != SynComponent::ISLOW);
    const bool use_i_slow = (syn_component != SynComponent::IFAST);

    typename ChemicalSynapseType::ConstructorArgs syn_args{};
    syn_args.params[ChemicalSynapseType::Esyn] = p["e_syn"].as<double>();

    if (use_i_fast)
    {
        syn_args.params[ChemicalSynapseType::gfast] = p["g_fast"].as<double>();
        syn_args.params[ChemicalSynapseType::sfast] = p["s_fast"].as<double>();
        syn_args.params[ChemicalSynapseType::Vfast] = p["v_fast"].as<double>();
    }
    else
    {
        syn_args.params[ChemicalSynapseType::gfast] = 0.0;
        syn_args.params[ChemicalSynapseType::sfast] = 0.0;
        syn_args.params[ChemicalSynapseType::Vfast] = 0.0;
    }

    if (use_i_slow)
    {
        syn_args.params[ChemicalSynapseType::gslow] = p["g_slow"].as<double>();
        syn_args.params[ChemicalSynapseType::k1] = p["k1"].as<double>();
        syn_args.params[ChemicalSynapseType::k2] = p["k2"].as<double>();
        syn_args.params[ChemicalSynapseType::sslow] = p["s_slow"].as<double>();
        syn_args.params[ChemicalSynapseType::Vslow] = p["v_slow"].as<double>();
    }
    else
    {
        syn_args.params[ChemicalSynapseType::gslow] = 0.0;
        syn_args.params[ChemicalSynapseType::k1] = 0.0;
        syn_args.params[ChemicalSynapseType::k2] = 0.0;
        syn_args.params[ChemicalSynapseType::sslow] = 0.0;
        syn_args.params[ChemicalSynapseType::Vslow] = 0.0;
    }

    syn_args.params[ChemicalSynapseType::mslow] = 0.0;

    ChemicalSynapseType synapse(create_hindmarsh_rose<Integrator>(true), NeuronType::x, model_neur, NeuronType::x, syn_args, 1);

    std::ofstream out(output_csv_path);
    if (syn_component == 0)
        out << "t,v_pre,v_post,i_fast\n";
    else if (syn_component == 1)
        out << "t,v_pre,v_post,i_slow\n";
    else
        out << "t,v_pre,v_post,i,i_fast,i_slow\n";

    const double dt = scaled_result->dt;
    const double *v_pre_sig_ptr = scaled_result->sig.data();
    const double *interpolated_points_ptr = scaled_result->interpolated_points.data();
    const size_t total_size = scaled_result->sig.size();
    const size_t points_factor = scaled_result->points_factor;

    constexpr auto i_enum = ChemicalSynapseType::i;
    constexpr auto i_fast_enum = ChemicalSynapseType::ifast;
    constexpr auto i_slow_enum = ChemicalSynapseType::islow;

    size_t interp_pts_counter = 0;
    size_t v_pre_sig_idx = 0;
    const size_t v_pre_sig_start_idx = static_cast<size_t>(stabilization_time / csv_step);

    // --- Fase de estabilización ---
    for (; v_pre_sig_idx < v_pre_sig_start_idx; v_pre_sig_idx++)
    {
        synapse.step(dt, v_pre_sig_ptr[v_pre_sig_idx], get_v_hindmarsh_rose(model_neur));
        const double i_val = synapse.get(i_enum);
        model_neur.add_synaptic_input(-i_val);
        model_neur.step(dt);

        for (size_t k = 1; k < points_factor; k++, interp_pts_counter++)
        {
            synapse.step(dt, interpolated_points_ptr[interp_pts_counter], get_v_hindmarsh_rose(model_neur));
            const double i_interp_val = synapse.get(i_enum);
            model_neur.add_synaptic_input(-i_interp_val);
            model_neur.step(dt);
        }
    }

    // --- Fase de simulación y registro ---
    for (; v_pre_sig_idx < total_size - 1; v_pre_sig_idx++)
    {
        const double v_pre = v_pre_sig_ptr[v_pre_sig_idx];

        synapse.step(dt, v_pre, get_v_hindmarsh_rose(model_neur));
        const double i_val = synapse.get(i_enum);

        model_neur.add_synaptic_input(-i_val);
        model_neur.step(dt);

        // Registro en el CSV
        if ((v_pre_sig_idx - v_pre_sig_start_idx) % subsampling == 0)
        {
            out << (v_pre_sig_idx - v_pre_sig_start_idx) * csv_step << "," << v_pre << "," << get_v_hindmarsh_rose(model_neur);
            if (syn_component == 0)
                out << "," << synapse.get(i_fast_enum) << "\n";
            else if (syn_component == 1)
                out << "," << synapse.get(i_slow_enum) << "\n";
            else
                out << "," << i_val << "," << synapse.get(i_fast_enum) << "," << synapse.get(i_slow_enum) << "\n";
        }

        // Sub-pasos de integración
        for (size_t k = 1; k < points_factor; k++, interp_pts_counter++)
        {
            synapse.step(dt, interpolated_points_ptr[interp_pts_counter], get_v_hindmarsh_rose(model_neur));
            const double i_interp_val = synapse.get(i_enum);
            model_neur.add_synaptic_input(-i_interp_val);
            model_neur.step(dt);
        }
    }

    // Último punto (sin sub-pasos posteriores)
    const double v_pre = v_pre_sig_ptr[v_pre_sig_idx];

    synapse.step(dt, v_pre, get_v_hindmarsh_rose(model_neur));
    const double i_val = synapse.get(i_enum);
    model_neur.add_synaptic_input(-i_val);
    model_neur.step(dt);

    if ((v_pre_sig_idx - v_pre_sig_start_idx) % subsampling == 0)
    {
        out << (v_pre_sig_idx - v_pre_sig_start_idx) * csv_step << "," << v_pre << "," << get_v_hindmarsh_rose(model_neur);
        if (syn_component == 0)
            out << "," << synapse.get(i_fast_enum) << "\n";
        else if (syn_component == 1)
            out << "," << synapse.get(i_slow_enum) << "\n";
        else
            out << "," << i_val << "," << synapse.get(i_fast_enum) << "," << synapse.get(i_slow_enum) << "\n";
    }

    out.close();
    return 0;
}
