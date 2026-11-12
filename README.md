[README.md](https://github.com/user-attachments/files/34729108/README.md)
# Implementación y evaluación de filtros digitales en Python

Este proyecto compara filtros FIR de ventana Hamming y filtros IIR Butterworth pasa bajos, pasa altos y pasa bandas. Incluye 18 configuraciones, gráficas y mediciones reproducibles.

## Ejecución

Requiere Python 3.12 o posterior. Desde esta carpeta:

```bash
python -m venv .venv
# Windows PowerShell:
.venv\Scripts\Activate.ps1
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
python filtros_digitales.py --salida resultados
python -m unittest -v
```

Si PowerShell impide activar el entorno, puede ejecutarse directamente `.venv\Scripts\python.exe` en lugar de `python`, sin cambiar la política del sistema.

## Señal de prueba

- Muestreo: 2000 Hz, 4 s, 8000 muestras, Nyquist: 1000 Hz.
- Tonos: 50 Hz (amplitud 1), 200 Hz (0.8), 400 Hz (0.6).
- Ruido blanco gaussiano: media 0, desviación estándar 0.5, semilla 2026.
- Cortes: pasa bajos 100 Hz, pasa altos 300 Hz, pasa bandas 150–250 Hz.
- Filtros principales: FIR de 201 coeficientes (orden 200) e IIR de orden de prototipo 4 (orden real 8 para pasa bandas).
- Variaciones: FIR de 51, 101 y 201 coeficientes; IIR de órdenes de prototipo 2, 4 y 6.

## Método de evaluación

Se aplica una sola pasada causal, con estado inicial nulo. No se utiliza `filtfilt`. Se excluyen los primeros 0.5 s de las métricas para evitar el arranque. Las gráficas temporales conservan la fase y el retardo.

La SNR de salida se calcula separando la salida del tono deseado de la salida de los otros tonos más ruido: 10 log10(P_deseada/P_no_deseada). Esto es posible porque la señal sintética y sus componentes son conocidas. La SNR no confunde el retardo de un filtro causal con ruido adicional. No se afirma que mida fidelidad respecto del tono original sin filtrar.

La respuesta se evalúa en 20001 frecuencias entre 0 y 1000 Hz. Criterio exploratorio: pérdida de paso <= 1 dB y rechazo mínimo >= 20 dB.

| Tipo | Banda de paso evaluada | Banda de rechazo evaluada |
|---|---|---|
| Pasa bajos | 0–70 Hz | 150–1000 Hz |
| Pasa altos | 400–1000 Hz | 0–200 Hz |
| Pasa bandas | 185–215 Hz | 0–100 Hz y 300–1000 Hz |

Los valores de corte no usan la misma convención: Butterworth tiene aproximadamente -3 dB en el corte; `firwin` usa aproximadamente -6 dB. La comparación principal mantiene los mismos cortes nominales y evalúa explícitamente bandas comunes. No compara diseños de igual costo ni de igual orden.

## Resultados principales

| Tipo | SNR entrada | SNR salida FIR | SNR salida IIR |
|---|---:|---:|---:|
| Pasa bajos | -1.76 dB | 12.94 dB | 12.60 dB |
| Pasa altos | -7.74 dB | 0.07 dB | -0.36 dB |
| Pasa bandas | -4.66 dB | 11.29 dB | 10.88 dB |

Los FIR de 201 coeficientes cumplen el criterio en las tres bandas y tienen retardo constante de 50 ms. Los IIR de prototipo 4 cumplen el criterio solo en pasa bandas. Los IIR de orden 6 sí lo cumplen en los tres casos. Un pasa altos deja pasar una fracción amplia del ruido blanco; reducir interferencia de baja frecuencia no elimina ruido situado en su banda útil.

El FIR pasa bandas de 51 coeficientes también cumple el criterio, con 12.5 ms de retardo y SNR 11.99 dB. Más coeficientes no garantizan mayor SNR para una señal específica: también importan el ancho de ruido efectivo y la forma de la respuesta.

## Publicación

El proyecto se entrega preparado para publicación manual. Consulta `SUBIR_A_GITHUB.md` para realizar `git commit` y `git push` al repositorio de destino.

## Archivos

- `filtros_digitales.py`: señal, diseño, filtrado, métricas y visualización.
- `test_filtros.py`: estabilidad, simetría FIR, convención de corte y amplitudes de tonos.
- `requirements.txt`: versiones utilizadas.
- `resultados/`: siete gráficas PNG, métricas CSV, coeficientes y resumen JSON.
- `Informe_filtros_digitales.pdf`: desarrollo y análisis.
- `Presentacion_filtros_digitales.pptx`: apoyo para exposición, sin portada.

## Guion para la exposición

1. Explicar qué tono interesa conservar en cada caso.
2. Mostrar las constantes de señal, `design`, `apply` y `measure`.
3. Ejecutar el script y abrir la respuesta en frecuencia y las comparaciones.
4. Comparar SNR, rechazo y retardo, y explicar el ruido que el pasa altos conserva.
5. Mostrar una configuración que no cumple y cómo cambia al aumentar el orden.

La sesión en línea y la exposición corresponden al estudiante. No se requiere MATLAB, Proteus ni Multisim.

## Referencias

SciPy. (s. f.). *butter*. https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html

SciPy. (s. f.). *firwin*. https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.firwin.html

SciPy. (s. f.). *sosfilt*. https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sosfilt.html

SciPy. (s. f.). *welch*. https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.welch.html
