#include "plugin-template.h"
#include <unistd.h>
#include <chrono>
#include <numeric>

extern "C" Plugin::Object* createRTXIPlugin(void) {
    return new PluginTemplate();
}

static DefaultGUIModel::variable_t vars[] = {
    {"ID Actualizacion", "Contador de ciclos NRT", DefaultGUIModel::STATE},
    {"Tamano Vector", "Tamano ajustado por NRT", DefaultGUIModel::STATE},
    {"Promedio Datos", "Calculado en RT", DefaultGUIModel::STATE},
    {"NRT Counter", "Contador incrementado en NRT", DefaultGUIModel::PARAMETER | DefaultGUIModel::UINTEGER},
};
static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

PluginTemplate::PluginTemplate(void)
    : DefaultGUIModel("Estructuras Dinamicas", ::vars, ::num_vars),
      idx_rt(0),
      gui_id_actualizacion(0),
      gui_tamano_vector(0),
      gui_promedio(0),
      gui_nrt_counter(0),
      hilo_activo(true),
      contador_nrt(0)
{
    // Inicialización básica de los vectores
    buffers[0].valores.assign(5, 0.0);
    buffers[1].valores.assign(5, 0.0);
    buffers[0].id_actualizacion = 0;
    buffers[1].id_actualizacion = 0;

    DefaultGUIModel::createGUI(vars, num_vars);
    customizeGUI();
    update(INIT);
    refresh();
}

PluginTemplate::~PluginTemplate(void) {
    hilo_activo = false; 
    if (hilo_nrt.joinable()) {
        hilo_nrt.join(); 
    }
}

void PluginTemplate::customizeGUI(void) {
    QGridLayout* layout = DefaultGUIModel::getLayout();
    QPushButton* btn = new QPushButton("Ejecutar Actualización NRT (10s)");
    layout->addWidget(btn, 0, 0); // Añadir botón a la rejilla superior
    QObject::connect(btn, SIGNAL(clicked()), this, SLOT(iniciarTareaNRT()));
}

void PluginTemplate::iniciarTareaNRT(void) {
    // Si hay un hilo previo que terminó, lo limpiamos antes de lanzar uno nuevo
    if (hilo_nrt.joinable()) {
        hilo_nrt.join();
    }
    hilo_nrt = std::thread(&PluginTemplate::tarea_pesada_nrt, this);
}

void PluginTemplate::execute(void) {
    int idx = idx_rt.load(std::memory_order_acquire);
    const std::vector<double>& v = buffers[idx].valores;
    
    if (!v.empty()) {
        double suma = std::accumulate(v.begin(), v.end(), 0.0);
        gui_promedio = suma / v.size();
        gui_tamano_vector = static_cast<double>(v.size());
        gui_id_actualizacion = static_cast<double>(buffers[idx].id_actualizacion);
    }
}

void PluginTemplate::update(DefaultGUIModel::update_flags_t flag) {
    switch (flag) {
        case INIT:
        case PERIOD:
            setState("ID Actualizacion", gui_id_actualizacion);
            setState("Tamano Vector", gui_tamano_vector);
            setState("Promedio Datos", gui_promedio);
            setParameter("NRT Counter", gui_nrt_counter);
            break;
        default:
            break;
    }
}

void PluginTemplate::tarea_pesada_nrt(void) {
    // Espera de 10 segundos interrumpible
    for (int i = 0; i < 100 && hilo_activo; ++i) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    
    if (!hilo_activo) return;

    contador_nrt++;
    gui_nrt_counter++;

    // Identificar el buffer inactivo
    int idx_escritura = (idx_rt.load() + 1) % 2;
    
    // --- Lógica de actualización única ---
    int nuevo_tamano = 5 + (contador_nrt % 10); 
    buffers[idx_escritura].valores.resize(nuevo_tamano);
    
    for (int j = 0; j < nuevo_tamano; ++j) {
        buffers[idx_escritura].valores[j] = (double)contador_nrt * (j + 1);
    }
    buffers[idx_escritura].id_actualizacion = contador_nrt;

    // Publicar datos al hilo RT
    idx_rt.store(idx_escritura, std::memory_order_release);
}