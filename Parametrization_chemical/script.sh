#!/usr/bin/env bash
# Script: Ejecuta generación de YAML, parametrización y plot
# Resultado final: la imagen tendrá el mismo nombre base que la salida .asc

set -euo pipefail
IFS=$'\n\t'

# Mover al directorio del script (espera ejecutarse desde la carpeta Parametrization_chemical)
cd "$(dirname "${BASH_SOURCE[0]}")"

echo "🔧 Starting full parametrization + plot pipeline"

# Activar entorno virtual si existe
if [ -f ".venv/bin/activate" ]; then
  # shellcheck source=/dev/null
  source .venv/bin/activate
  echo "✅ Virtualenv activated (.venv)"
else
  echo "⚠️  No .venv found — using system Python (asegúrate del entorno correcto)"
fi

# 1) Limpiar yaml_configs
echo "
--- Cleaning yaml_configs/ ---"
rm -rf yaml_configs/*

# 2) Generar YAML
echo "\n--- Running yaml_generator.py ---"
python yaml_generator.py

# Detectar YAML generado (el más reciente si hay varios)
yaml_file=$(ls -1t yaml_configs/*.yaml 2>/dev/null | head -n1 || true)
if [ -z "${yaml_file}" ]; then
  echo "❌ No YAML generated. Abort."
  exit 1
fi
base_name=$(basename "${yaml_file}" .yaml)
asc_file="yaml_configs/${base_name}.asc"

echo "Saved YAML: ${yaml_file}"

# 3) Ejecutar parametrización (procesa todos los YAMLs en yaml_configs)
echo "\n--- Running run_parametrization.py ---"
python run_parametrization.py

# 4) Verificar salida .asc (buscar por base_name si no existe exactamente)
if [ ! -f "${asc_file}" ]; then
  asc_file_found=$(ls -1t yaml_configs/*.asc 2>/dev/null | grep "${base_name}" | head -n1 || true)
  if [ -n "${asc_file_found}" ]; then
    asc_file="${asc_file_found}"
  else
    echo "❌ No .asc output found for ${base_name}. Contents of yaml_configs/:"
    ls -l yaml_configs || true
    exit 1
  fi
fi

echo "Output .asc: ${asc_file}"

# 5) Plotear (plot.py genera images/<base>.eps)
echo "\n--- Running plot.py on ${asc_file} ---"
python plot.py "${asc_file}" 1

image_file="images/${base_name}.eps"
if [ -f "${image_file}" ]; then
  echo "✅ Plot creado: ${image_file}"
else
  echo "⚠️ plot.py no generó ${image_file} — revisa la salida de plot.py"
fi

echo "\n🎉 Pipeline completado. Archivo final de imagen: ${image_file}"
exit 0
