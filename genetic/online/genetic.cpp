#include "bidirectional_chemical_synapse_genetic.h"
#include <chrono>
#include <iostream>

void BidirectionalChemicalSynapseGenetic::NRT_genetic(double period_t, double dt_factor_t, double v_max_1_t, double v_min_1_t, double v_max_2_t, double v_min_2_t)
{
    const bool use_syn_21 = use_i_fast_21 || use_i_slow_21;
    const bool use_syn_12 = use_i_fast_12 || use_i_slow_12;

    if (use_syn_21 || use_syn_12)
    {
        const size_t curr_synapse_idx = synapse_idx.load(std::memory_order_relaxed);
        const size_t new_synapse_idx = 1 - curr_synapse_idx;

        params_12[new_synapse_idx].e_syn += params_12[curr_synapse_idx].e_syn + 0.1;
        params_12[new_synapse_idx].g_fast += params_12[curr_synapse_idx].g_fast + 0.1;
        params_12[new_synapse_idx].s_fast += params_12[curr_synapse_idx].s_fast + 0.1;
        params_12[new_synapse_idx].v_fast += params_12[curr_synapse_idx].v_fast + 0.1;
        params_12[new_synapse_idx].g_slow += params_12[curr_synapse_idx].g_slow + 0.1;
        params_12[new_synapse_idx].k1 += params_12[curr_synapse_idx].k1 + 0.1;
        params_12[new_synapse_idx].k2 += params_12[curr_synapse_idx].k2 + 0.1;
        params_12[new_synapse_idx].s_slow += params_12[curr_synapse_idx].s_slow + 0.1;
        params_12[new_synapse_idx].v_slow += params_12[curr_synapse_idx].v_slow + 0.1;

        params_21[new_synapse_idx].e_syn += params_21[curr_synapse_idx].e_syn + 0.1;
        params_21[new_synapse_idx].g_fast += params_21[curr_synapse_idx].g_fast + 0.1;
        params_21[new_synapse_idx].s_fast += params_21[curr_synapse_idx].s_fast + 0.1;
        params_21[new_synapse_idx].v_fast += params_21[curr_synapse_idx].v_fast + 0.1;
        params_21[new_synapse_idx].g_slow += params_21[curr_synapse_idx].g_slow + 0.1;
        params_21[new_synapse_idx].k1 += params_21[curr_synapse_idx].k1 + 0.1;
        params_21[new_synapse_idx].k2 += params_21[curr_synapse_idx].k2 + 0.1;
        params_21[new_synapse_idx].s_slow += params_21[curr_synapse_idx].s_slow + 0.1;
        params_21[new_synapse_idx].v_slow += params_21[curr_synapse_idx].v_slow + 0.1;

        synapse_idx.store(new_synapse_idx, std::memory_order_release);

        QMetaObject::invokeMethod(this, "update_params_gui", Qt::QueuedConnection);

        if (stabilization_time > 0.0)
        {
            std::this_thread::sleep_for(std::chrono::duration<double>(stabilization_time));
        }

        num_elements = (int)(evaluation_time / (period_t * 1e-3));

        if (use_syn_12)
        {
            v_sig_1.resize(num_elements);
            if (use_i_fast_12)
                i_fast_sig_12.resize(num_elements);
            if (use_i_slow_12)
                i_slow_sig_12.resize(num_elements);
        }

        if (use_syn_21)
        {
            v_sig_2.resize(num_elements);
            if (use_i_fast_21)
                i_fast_sig_21.resize(num_elements);
            if (use_i_slow_21)
                i_slow_sig_21.resize(num_elements);
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
            print_first("v_sig_1", v_sig_1);
            if (use_i_fast_12)
                print_first("i_fast_sig_12", i_fast_sig_12);
            if (use_i_slow_12)
                print_first("i_slow_sig_12", i_slow_sig_12);
        }

        if (use_syn_21)
        {
            print_first("v_sig_2", v_sig_2);
            if (use_i_fast_21)
                print_first("i_fast_sig_21", i_fast_sig_21);
            if (use_i_slow_21)
                print_first("i_slow_sig_21", i_slow_sig_21);
        }
    }

    QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
}
