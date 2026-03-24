#include "bidirectional_chemical_synapse_genetic.h"
#include <chrono>

void BidirectionalChemicalSynapseGenetic::NRT_genetic(void)
{
    bool use_syn_21 = use_i_fast_21 || use_i_slow_21;
    bool use_syn_12 = use_i_fast_12 || use_i_slow_12;
    bool use_syn = use_syn_21 || use_syn_12;

    if (use_syn)
    {
        synapse_lock.acquire();
        m_slow_12 = 0.0;
        m_slow_21 = 0.0;
        for (int i = 0; i < SP_COUNT; ++i)
        {
            params_12[i] += 1.0;
            params_21[i] += 1.0;
        }
        synapse_lock.release();

        if (stabilization_time > 0.0)
        {
            std::this_thread::sleep_for(std::chrono::duration<double>(stabilization_time));
        }

        num_elements = (int)(evaluation_time * freq);

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

        QMetaObject::invokeMethod(this, "update_params_gui", Qt::QueuedConnection);
    }

    QMetaObject::invokeMethod(this, "stop_genetic_event_async", Qt::QueuedConnection);
}
