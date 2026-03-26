#include "bidirectional_chemical_synapse_genetic.h"
#include <chrono>
#include <iostream>

void BidirectionalChemicalSynapseGenetic::NRT_genetic(double thread_freq, double scale_21_t, double offset_21_t, double scale_12_t, double offset_12_t, unsigned int i_ranges_from_neuron_t, double i_max_21_t, double i_min_21_t, double i_max_12_t, double i_min_12_t, double max_1_t, double min_1_t)
{
    const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;
    const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;

    if (use_syn_21 || use_syn_12)
    {
        const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_relaxed);
        const size_t new_synapse_idx = 1 - curr_synapse_idx;
        m_slow_12[new_synapse_idx] = 0.0;
        m_slow_21[new_synapse_idx] = 0.0;

        params_12[new_synapse_idx].e_syn += 0.1;
        params_12[new_synapse_idx].g_fast += 0.1;
        params_12[new_synapse_idx].s_fast += 0.1;
        params_12[new_synapse_idx].v_fast += 0.1;
        params_12[new_synapse_idx].g_slow += 0.1;
        params_12[new_synapse_idx].k1 += 0.1;
        params_12[new_synapse_idx].k2 += 0.1;
        params_12[new_synapse_idx].s_slow += 0.1;
        params_12[new_synapse_idx].v_slow += 0.1;

        params_21[new_synapse_idx].e_syn += 0.1;
        params_21[new_synapse_idx].g_fast += 0.1;
        params_21[new_synapse_idx].s_fast += 0.1;
        params_21[new_synapse_idx].v_fast += 0.1;
        params_21[new_synapse_idx].g_slow += 0.1;
        params_21[new_synapse_idx].k1 += 0.1;
        params_21[new_synapse_idx].k2 += 0.1;
        params_21[new_synapse_idx].s_slow += 0.1;
        params_21[new_synapse_idx].v_slow += 0.1;
        synapse_idx.store(new_synapse_idx, std::memory_order_release);

        QMetaObject::invokeMethod(this, "update_params_gui", Qt::QueuedConnection);

        if (stabilization_time > 0.0)
        {
            std::this_thread::sleep_for(std::chrono::duration<double>(stabilization_time));
        }

        num_elements = (int)(evaluation_time * thread_freq);

        if (use_syn_12)
        {
            v1_scaled_sig.resize(num_elements);
            v2_sig.resize(num_elements);
            if (use_i_fast_12)
                i_fast_sig_12.resize(num_elements);
            if (use_i_slow_12)
                i_slow_12.resize(num_elements);
        }

        if (use_syn_21)
        {
            v2_scaled_sig.resize(num_elements);
            v1_sig.resize(num_elements);
            if (use_i_fast_21)
                i_fast_sig_21.resize(num_elements);
            if (use_i_slow_21)
                i_slow_21.resize(num_elements);
        }

        storing_idx = 0;
        RT_storing.store(true, std::memory_order_release);

        while (RT_storing.load(std::memory_order_acquire))
        {
            if (stop_genetic.load(std::memory_order_relaxed))
            {
                RT_storing.store(false, std::memory_order_relaxed);
                QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
                return;
            }
            std::this_thread::sleep_for(std::chrono::duration<double>(0.01));
        }

        QMetaObject::invokeMethod(this, "set_generations_completed", Qt::QueuedConnection, Q_ARG(double, generations_completed + 1));
        QMetaObject::invokeMethod(this, "set_individuals_completed", Qt::QueuedConnection, Q_ARG(double, individuals_completed + 1));

        const size_t act_scaling_factors_21_idx = scaling_factors_21_idx.load(std::memory_order_acquire);
        std::cout << "Final scale 21: " << scale_21[act_scaling_factors_21_idx] << ", offset 21: " << offset_21[act_scaling_factors_21_idx] << "\n";

        const size_t act_scaling_factors_12_idx = scaling_factors_12_idx.load(std::memory_order_acquire);
        std::cout << "Final scale 12: " << scale_12[act_scaling_factors_12_idx] << ", offset 12: " << offset_12[act_scaling_factors_12_idx] << "\n";

        auto print_first = [&](const char *name, const univector<double> &v)
        {
            int n = (int)std::min<size_t>(v.size(), 10);
            std::cout << name << " [" << v.size() << "] -> ";
            for (int i = 0; i < n; i++)
            {
                std::cout << v[i];
                if (i + 1 < n)
                    std::cout << ", ";
            }
            std::cout << "\n";
        };

        if (use_syn_12)
        {
            print_first("v1_scaled_sig", v1_scaled_sig);
            print_first("v2_sig", v2_sig);
            if (use_i_fast_12)
                print_first("i_fast_sig_12", i_fast_sig_12);
            if (use_i_slow_12)
                print_first("i_slow_12", i_slow_12);
        }

        if (use_syn_21)
        {
            print_first("v2_scaled_sig", v2_scaled_sig);
            print_first("v1_sig", v1_sig);
            if (use_i_fast_21)
                print_first("i_fast_sig_21", i_fast_sig_21);
            if (use_i_slow_21)
                print_first("i_slow_21", i_slow_21);
        }
    }

    QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
}
