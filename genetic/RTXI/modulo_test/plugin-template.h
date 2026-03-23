#include <default_gui_model.h>
#include <thread>
#include <atomic>
#include <vector>
#include <QPushButton>
#include <QGridLayout>

// Estructura que contiene el vector dinámico y un ID de versión
struct DatosComplejos {
    std::vector<double> valores;
    int id_actualizacion;
};

class PluginTemplate : public DefaultGUIModel
{
    Q_OBJECT

public:
    PluginTemplate(void);
    virtual ~PluginTemplate(void);
    void execute(void); // Hilo de Tiempo Real

protected:
    virtual void update(DefaultGUIModel::update_flags_t); // Hilo de Interfaz

private slots:
    // Slot para iniciar el hilo desde la GUI
    void iniciarTareaNRT(void);

private:
    // Mecanismo de comunicación Lock-Free
    DatosComplejos buffers[2];
    std::atomic<int> idx_rt; // Índice que el hilo RT lee actualmente

    // Variables persistentes para la GUI
    double gui_id_actualizacion;
    double gui_tamano_vector;
    double gui_promedio;

    // Nuevo parámetro NRT (incrementado desde el hilo no real-time)
    int gui_nrt_counter;

    // Control del hilo No-Real-Time
    std::atomic<bool> hilo_activo;
    std::thread hilo_nrt;
    int contador_nrt; // Contador persistente para los cálculos

    void tarea_pesada_nrt(void);
    void customizeGUI(void);
};