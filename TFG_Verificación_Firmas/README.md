# TFG: Verificación de Firmas en Iniciativas Legislativas Populares (Castilla y León)

## 📌 Descripción
Este repositorio contiene el código desarrollado para el Trabajo de Fin de Grado en Estadística.  
El objetivo es diseñar metodologías estadísticas y de Machine Learning para optimizar la verificación de firmas en ILPs, aplicado al caso de Castilla y León.

## 📂 Estructura del repositorio
- `codigo/python/`: Notebooks de simulación de datos y modelos de Machine Learning (supervisados y no supervisados).
- `codigo/matlab/`: Scripts de diseños muestrales (estratificado, conglomerados, secuencial de Wald) y planes de aceptación.
- `datos/`: Datasets simulados (10% y 20% de inválidas) generados por los notebooks.
- `resultados/`: Figuras y tablas generadas durante el estudio.

## ⚙️ Requisitos e instalación
- **Python 3.8+**: Instalar dependencias con `pip install -r requirements.txt`
- **MATLAB R2020a+**: Ejecutar los scripts `.m` directamente (requiere Statistics and Machine Learning Toolbox).

## 🚀 Cómo ejecutar el código
1. **Machine Learning**: Navega a `codigo/python/` y ejecuta en orden:
   - `GeneraciónDeDFYValidación.ipynb` (genera los CSV).
   - `TAA_TFG.ipynb` (entrena y evalúa los modelos).
2. **Muestreo**: Abre MATLAB, ve a `codigo/matlab/` y ejecuta `muestreo_TFG.m`.

## 👨‍💻 Autores
- Álvaro Gómez Jorge
- Tutor: Jesús Alberto Tapia García

## 📄 Licencia
Este proyecto está bajo la licencia MIT.