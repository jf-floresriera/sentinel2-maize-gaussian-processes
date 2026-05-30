<div align="center">

# 🌿 GPR Phenology — Finca La Esperanza
## Meta, Colombia · 2024

![Python](https://img.shields.io/badge/Python-3.11-3776AB?style=for-the-badge&logo=python&logoColor=white)
![GEE](https://img.shields.io/badge/Google_Earth_Engine-4285F4?style=for-the-badge&logo=google&logoColor=white)
![Sentinel-2](https://img.shields.io/badge/Sentinel--2-SR_Harmonized-00b4d8?style=for-the-badge)
![Tile](https://img.shields.io/badge/Tile-18NZK-ff6b35?style=for-the-badge)
![Año](https://img.shields.io/badge/Año-2024-success?style=for-the-badge)
![License](https://img.shields.io/badge/License-MIT-brightgreen?style=for-the-badge)

**Monitoreo de fenología de cultivos mediante modelos de**
**Regresión por Procesos Gaussianos (GPR) sobre imágenes Sentinel-2**
**procesadas en Google Earth Engine**

| 📍 Sitio | 🛰️ Satélite | 🗓️ Período | 🌱 Variables | 📐 Resolución |
|:---:|:---:|:---:|:---:|:---:|
| Finca La Esperanza, Meta · Colombia | Sentinel-2 SR Harmonized | Ene – Dic 2024 | LAI · FVC · laiCab | 20 m/px |

</div>

---

## 🗺️ Área de Interés

El área de estudio es la **Finca La Esperanza**, Meta, Colombia (~4.33° N, 72.01° O).
El AOI excluye el lote de la vivienda: `Area_Total.difference(Lote_Casa)` → **~150 ha netas de cultivo**.

> Máscara de nubes con SCL (clases 1–3 y 6–11) + QA60 (bits 10–11).
> `.clip(AOI)` aplicado **antes** del modelo GPR para eficiencia computacional.

---

## 🔬 Variables de Vegetación

| Variable | Nombre Completo | Unidades | Rango Típico | Resolución |
|:---:|:---|:---:|:---:|:---:|
| **LAI** | Leaf Area Index | m²/m² | 0 – 7 | 20 m |
| **FVC** | Fractional Vegetation Cover | fracción | 0 – 1 | 20 m |
| **laiCab** | Canopy Chlorophyll Content | µg/cm² | 0 – 600 | 20 m |

---

## 🌱 Métricas Fenológicas (LSP)

Calculadas con **doble logística** sobre la serie gap-filled, umbral 30% amplitud relativa:

| Banda `.tif` | Métrica | Descripción |
|:---:|:---:|:---|
| `sos` | **S**tart **o**f **S**eason | Inicio del ciclo de crecimiento (DOY) |
| `eos` | **E**nd **o**f **S**eason | Fin del ciclo de crecimiento (DOY) |
| `pos` | **P**eak **o**f **S**eason | Máximo índice de vegetación (DOY) |
| `los` | **L**ength **o**f **S**eason | Duración del ciclo activo (días) |

---

## 🔄 Pipeline de Procesamiento

```
🛰️  Sentinel-2 SR Harmonized  (COPERNICUS/S2_SR_HARMONIZED)
           │   Tile 18NZK · Bandas B2–B12 · 2024-01-01 → 2024-12-31
           │   Máscara SCL+QA60 · .clip(AOI) · CLOUDY_PIXEL_PERCENTAGE < 80
           ▼
📡  PASO 1 — GPR Predicted Mean               1_Paso_GEE_Adaptado.js
           │   Modelo GPR (LAI · FVC · laiCab) · 20 m · EPSG:4326
           ▼
      📁  GPR_LAI_2024/  ·  GPR_FVC_2024/  ·  GPR_laiCab_2024/
           ▼
🔄  PASO 2 — GPR Gap-filling ±30 días         2_Paso_GEE_Adaptado.js
           │   Kernel GPR · ell2_ts · sigf_ts · sign_ts
           ▼
      📁  GF_LAI_2024/   ·  GF_FVC_2024/   ·  GF_laiCab_2024/
           ▼
🌿  PASO 3 — LSP Generation                   3_Paso_GEE_Adaptado.js
           │   Doble Logística · SOS · EOS · POS · LOS
           ▼
      🗺️  LSP_LAI_2024.tif · LSP_FVC_2024.tif · LSP_laiCab_2024.tif
           ▼
🐍  ANÁLISIS LOCAL EN PYTHON
    a_descargar_GEE.py  ·  b_calidad.py  ·  c_validacion.py  ·  d_graficas_articulo5.py
```

---

## 📊 Resultados — Control de Calidad

> Carpeta: `figuras_analysis/0_calidad/`

### Cobertura por variable — Gap-filled

| LAI | FVC | laiCab |
|:---:|:---:|:---:|
| ![](figuras_analysis/0_calidad/calidad_cobertura_GF_LAI.png) | ![](figuras_analysis/0_calidad/calidad_cobertura_GF_FVC.png) | ![](figuras_analysis/0_calidad/calidad_cobertura_GF_laiCab.png) |

### Cobertura por variable — GPR (sin gap-filling)

| LAI | FVC | laiCab |
|:---:|:---:|:---:|
| ![](figuras_analysis/0_calidad/calidad_cobertura_GPR_LAI.png) | ![](figuras_analysis/0_calidad/calidad_cobertura_GPR_FVC.png) | ![](figuras_analysis/0_calidad/calidad_cobertura_GPR_laiCab.png) |

### Resumen boxplot de calidad

<div align="center">

![Boxplot calidad](figuras_analysis/0_calidad/calidad_resumen_boxplot.png)

</div>

---

## 📈 Resultados — Validación

> Carpeta: `figuras_analysis/2_validacion_v3/`

### Series temporales con percentiles

| LAI | FVC | laiCab |
|:---:|:---:|:---:|
| ![](figuras_analysis/2_validacion_v3/temporal_percentiles_LAI.png) | ![](figuras_analysis/2_validacion_v3/temporal_percentiles_FVC.png) | ![](figuras_analysis/2_validacion_v3/temporal_percentiles_laiCab.png) |

### Métricas temporales

| LAI | FVC | laiCab |
|:---:|:---:|:---:|
| ![](figuras_analysis/2_validacion_v3/metricas_temporal_LAI.png) | ![](figuras_analysis/2_validacion_v3/metricas_temporal_FVC.png) | ![](figuras_analysis/2_validacion_v3/metricas_temporal_laiCab.png) |

### Scatter plots balanceados

| LAI | FVC | laiCab |
|:---:|:---:|:---:|
| ![](figuras_analysis/2_validacion_v3/scatter_balanced_LAI.png) | ![](figuras_analysis/2_validacion_v3/scatter_balanced_FVC.png) | ![](figuras_analysis/2_validacion_v3/scatter_balanced_laiCab.png) |

### Scatter + Mapas combinados

<div align="center">

![Scatter maps](figuras_analysis/2_validacion_v3/FigA_scatter_maps.png)

</div>

---

## 🖼️ Figuras Finales para Publicación

> Carpeta: `figuras_analysis/3_final_figure/`

<div align="center">

| | |
|:---:|:---:|
| ![Figura 1](figuras_analysis/3_final_figure/figura1.png) | ![Figura 2](figuras_analysis/3_final_figure/figura2.png) |
| ![Figura 3](figuras_analysis/3_final_figure/figura3.png) | ![Figura 4](figuras_analysis/3_final_figure/figura4.png) |
| ![Figura 5](figuras_analysis/3_final_figure/figura5.png) | ![Figura 6](figuras_analysis/3_final_figure/figura6.png) |
| ![Figura 7](figuras_analysis/3_final_figure/figura7.png) | ![Figura 8](figuras_analysis/3_final_figure/figura8.png) |

</div>

---

## 📁 Estructura del Repositorio

```
GPR-Phenology-FincaLaEsperanza/
│
├── 📂 GEE_Downloads_tiff/
│   ├── 📂 GPR_LAI_2024/         · LAI por fecha, sin gap-filling    (Paso 1)
│   ├── 📂 GPR_FVC_2024/         · FVC por fecha, sin gap-filling    (Paso 1)
│   ├── 📂 GPR_laiCab_2024/      · laiCab por fecha, sin gap-filling (Paso 1)
│   ├── 📂 GF_LAI_2024/          · LAI gap-filled por fecha          (Paso 2)
│   ├── 📂 GF_FVC_2024/          · FVC gap-filled por fecha          (Paso 2)
│   ├── 📂 GF_laiCab_2024/       · laiCab gap-filled por fecha       (Paso 2)
│   ├── 🗺️ LSP_LAI_2024.tif      · SOS/EOS/POS/LOS · LAI            (Paso 3)
│   ├── 🗺️ LSP_FVC_2024.tif      · SOS/EOS/POS/LOS · FVC            (Paso 3)
│   └── 🗺️ LSP_laiCab_2024.tif   · SOS/EOS/POS/LOS · laiCab         (Paso 3)
│
├── 📂 code/
│   ├── 📂 CODIGO_GEE_original_msalinero/
│   ├── 🌍 1_Paso_GEE_Adaptado.js
│   ├── 🌍 2_Paso_GEE_Adaptado.js
│   ├── 🌍 3_Paso_GEE_Adaptado.js
│   ├── 🐍 a_descargar_GEE.py
│   ├── 🐍 b_calidad.py
│   ├── 🐍 c_validacion.py
│   └── 🐍 d_graficas_articulo5.py
│
├── 📂 figuras_analysis/
│   ├── 📂 0_calidad/
│   │   ├── calidad_cobertura_GF_FVC.png
│   │   ├── calidad_cobertura_GF_LAI.png
│   │   ├── calidad_cobertura_GF_laiCab.png
│   │   ├── calidad_cobertura_GPR_FVC.png
│   │   ├── calidad_cobertura_GPR_LAI.png
│   │   ├── calidad_cobertura_GPR_laiCab.png
│   │   ├── calidad_resumen_boxplot.png
│   │   └── resumen_calidad_imagenes.csv
│   ├── 📂 1_validacion/
│   ├── 📂 2_validacion_v3/
│   │   ├── FigA_scatter_maps.png
│   │   ├── metricas_temporal_FVC.png / LAI / laiCab
│   │   ├── scatter_balanced_FVC.png / LAI / laiCab
│   │   ├── temporal_percentiles_FVC.png / LAI / laiCab
│   │   └── *.csv  (metricas, puntos_100, resumen_global)
│   └── 📂 3_final_figure/
│       ├── figura1.png  …  figura8.png
│
├── 📄 CONTEXT.md
├── 📄 README.md
├── 📄 environment.yml
└── 📄 LICENSE
```

---

## ⚙️ Instalación

```bash
# 1. Clonar
git clone https://github.com/<tu-usuario>/GPR-Phenology-FincaLaEsperanza.git
cd GPR-Phenology-FincaLaEsperanza

# 2. Entorno conda
conda env create -f environment.yml
conda activate gpr-phenology

# 3. Autenticar GEE
earthengine authenticate
earthengine set_project wide-origin-466923-d8

# 4. Aceptar repo de modelos GPR en GEE
# https://code.earthengine.google.com/?accept_repo=users/msalinero85/GPRPhenologyDemos
```

---

## 🚀 Ejecución

### Scripts GEE (en orden — Code Editor)

```
① Crear 6 ImageCollections en GEE Assets:
   GPR_LAI_2024  GPR_FVC_2024  GPR_laiCab_2024
   GF_LAI_2024   GF_FVC_2024   GF_laiCab_2024

② code/1_Paso_GEE_Adaptado.js  →  Tasks → Run All  (~50–80 tasks/variable)
③ Esperar Completed → code/2_Paso_GEE_Adaptado.js  →  Tasks → Run All
④ Esperar Completed → code/3_Paso_GEE_Adaptado.js  →  Tasks → Run All
```

### Scripts Python (en orden)

```bash
python code/a_descargar_GEE.py        # Descarga assets → GEE_Downloads_tiff/
python code/b_calidad.py              # Control de calidad
python code/c_validacion.py           # Validación
python code/d_graficas_articulo5.py   # Figuras → figuras_analysis/3_final_figure/
```

---

## 🔧 Configuración GEE

| Parámetro | Valor |
|:---|:---|
| Proyecto | `wide-origin-466923-d8` |
| Colección S2 | `COPERNICUS/S2_SR_HARMONIZED` |
| Tile | `18NZK` |
| Período | `2024-01-01` → `2024-12-31` |
| Escala | 20 m/px · `EPSG:4326` |
| Filtro nubes | `CLOUDY_PIXEL_PERCENTAGE < 80` |
| Máscara SCL | Clases 1, 2, 3, 6, 7, 8, 9, 10, 11 |
| Máscara QA60 | Bits 10 y 11 |
| Ventana gap-filling | ±30 días |
| Umbral LSP | 0.3 (30% amplitud relativa) |
| maxPixels Paso 1 | `1e9` |
| maxPixels Paso 2–3 | `1e13` |

---

## 📦 .gitignore recomendado

```gitignore
GEE_Downloads_tiff/
__pycache__/
*.pyc
*.pyo
.ipynb_checkpoints/
.env
*.log
```

---

<div align="center">
<sub>🛰️ Teledetección Agrícola · Google Earth Engine · GPR · Meta, Colombia · 2024–2026</sub>
</div>
