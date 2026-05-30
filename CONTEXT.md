# CONTEXT.md — Proyecto GPR Phenology Finca La Esperanza

> Documento de contexto técnico para transferencia a otro agente IA o colaborador.
> Generado: 2026-05-29

---

## 1. Investigador

| Campo | Detalle |
|---|---|
| **Nombre** | Jesus Enrique Flores Riera |
| **Institución** | Universidad Nacional de Colombia |
| **Ubicación** | Bogotá, Colombia |
| **Perfil** | Investigador en teledetección agrícola, GEE, Python, ML (YOLO, GPR) |

---

## 2. Referencia metodológica base

| Campo | Detalle |
|---|---|
| **Artículo** | "Monitoring Cropland Phenology on Google Earth Engine Using Gaussian Process Regression" |
| **Autores** | Salinero-Delgado et al. (2022) |
| **Revista** | Remote Sensing, MDPI, Vol. 14, No. 1, Art. 146 |
| **DOI** | https://doi.org/10.3390/rs14010146 |
| **Repo GEE original** | `users/msalinero85/GPRPhenologyDemos` |
| **URL aceptar repo** | https://code.earthengine.google.com/?accept_repo=users/msalinero85/GPRPhenologyDemos |

### Scripts del repositorio original (msalinero85)

| Script | Función |
|---|---|
| `S2BOAModels` | Hiperparámetros GPR (LAI, FVC, laiCab) exportados desde ARTMO/MATLAB a GEE |
| `visualization` | Paletas de color y rangos min/max por variable |
| `PhenologyFunctions` | Doble logística + `addDoy` + `get_Double_Logistic_Params` para métricas LSP |
| `GPRPredictedMean` | Demo: mapa GPR media predicha para tile 30TUM completo |
| `GPRGapfilling` | Demo: gap-filling temporal GPR sobre tile completo |
| `LSPGeneration` | Demo: métricas SOS/EOS/POS/LOS desde colección gap-filled |

---

## 3. Área de estudio

- **Sitio:** Finca La Esperanza, Meta, Colombia
- **Tile Sentinel-2:** `18NZK`
- **Período:** año completo 2024 (`2024-01-01` → `2024-12-31`)
- **AOI:** `Area_Total.difference(Lote_Casa)` — finca completa excluyendo el área de la vivienda
- **Proyecto GEE:** `wide-origin-466923-d8`
- **Asset root GEE:** `projects/wide-origin-466923-d8/assets/`

### Coordenadas AOI

```javascript
var Area_Total = ee.Geometry.Polygon([[
  [-72.02911590018898, 4.323774005392068],
  [-72.02320164960223, 4.320341858927772],
  [-71.99426061433154, 4.337725140154523],
  [-71.99372417252856, 4.338110271794471],
  [-71.99303752702075, 4.338666572705499],
  [-71.99149257462817, 4.340207096166426],
  [-71.99072009843188, 4.340870376132908],
  [-71.98967671912509, 4.341426675010239],
  [-72.00730581274610, 4.349063111509516],
  [-72.00862576723136, 4.353153345576473],
  [-72.01213946104087, 4.356828060080524],
  [-72.01526010903413, 4.352303621556992],
  [-72.01894546422060, 4.346895815199964],
  [-72.02098296949191, 4.341000332940370],
  [-72.02443983146695, 4.334085615058952],
  [-72.02447201797513, 4.333775368489792],
  [-72.02456857749966, 4.333529310775441],
  [-72.02911590018898, 4.323774005392068]
]]);

var Lote_Casa = ee.Geometry.Polygon([[
  [-72.02527036942277, 4.327000978587103],
  [-72.02123566610653, 4.329945886052444],
  [-72.02102108938534, 4.329988678906903],
  [-72.02003403646786, 4.329753318177450],
  [-72.01848908407528, 4.329603543129691],
  [-72.01829596502621, 4.329582146691876],
  [-72.01656862242062, 4.328672797525615],
  [-72.01685830099423, 4.326244294986618],
  [-72.01699704976303, 4.325322370602150],
  [-72.01762468667252, 4.324659076998467],
  [-72.01774806828720, 4.324616283842795],
  [-72.01786608548386, 4.324659076998467],
  [-72.01965780110581, 4.325525637880559],
  [-72.01980264039261, 4.325616573224313],
  [-72.02000648827774, 4.325718206830860],
  [-72.02018351407273, 4.325793094742748],
  [-72.02030689568741, 4.325916124867642],
  [-72.02041954846604, 4.326097995450447],
  [-72.02046782822830, 4.326124741120704],
  [-72.02074141354782, 4.326172883324792],
  [-72.02128858418686, 4.326199628992403],
  [-72.02208251805527, 4.325456099080927],
  [-72.02527036942277, 4.327000978587103]
]]);

var AOI    = Area_Total.difference({ right: Lote_Casa, maxError: 1 });
var parcel = AOI.bounds();
```

---

## 4. Pipeline implementado (3 pasos)

```
Sentinel-2 COPERNICUS/S2_SR_HARMONIZED (tile 18NZK, 2024)
        │
        ▼
  PASO 1 — GPR Predicted Mean  [1_Paso_GEE_Adaptado.js]
  maskS2cloud_and_water → .divide(scaleFactor) → .clip(parcel) → makeGPR(model)
  Export por imagen → assets GPR_LAI_2024/, GPR_FVC_2024/, GPR_laiCab_2024/
  Descargado como → GEE_Downloads_tiff/GPR_LAI_2024/, GPR_FVC_2024/, GPR_laiCab_2024/
        │
        ▼
  PASO 2 — GPR Gap-filling temporal ±30 días  [2_Paso_GEE_Adaptado.js]
  makeGapFill() con hiperparámetros ell2_ts, sigf_ts, sign_ts
  Export por imagen → assets GF_LAI_2024/, GF_FVC_2024/, GF_laiCab_2024/
  Descargado como → GEE_Downloads_tiff/GF_LAI_2024/, GF_FVC_2024/, GF_laiCab_2024/
        │
        ▼
  PASO 3 — LSP Generation (doble logística)  [3_Paso_GEE_Adaptado.js]
  get_Double_Logistic_Params(gapfilledCol, 'gapfilled', 0.3)
  Métricas: SOS, EOS, POS, LOS (día del año)
  Descargado como → GEE_Downloads_tiff/LSP_LAI_2024.tif, LSP_FVC_2024.tif, LSP_laiCab_2024.tif
        │
        ▼
  ANÁLISIS LOCAL (Python)
  a_descargar_GEE.py   → descarga assets de GEE a disco con geemap
  b_calidad.py         → control de calidad de los GeoTIFFs
  c_validacion.py      → validación de resultados
  d_graficas_articulo5.py → figuras para publicación
```

---

## 5. Variables procesadas

| Variable | Descripción | Unidades |
|---|---|---|
| `LAI` | Leaf Area Index | m²/m² |
| `FVC` | Fractional Vegetation Cover | fracción 0–1 |
| `laiCab` | Canopy Chlorophyll Content | µg/cm² |

## 6. Métricas LSP generadas (bandas de los .tif LSP)

| Banda | Descripción | Unidades |
|---|---|---|
| `sos` | Start of Season | día del año (DOY) |
| `eos` | End of Season | día del año (DOY) |
| `pos` | Peak of Season | día del año (DOY) |
| `los` | Length of Season | días |

---

## 7. Configuración GEE

| Parámetro | Valor |
|---|---|
| Colección S2 | `COPERNICUS/S2_SR_HARMONIZED` |
| Tile | `18NZK` |
| Escala export | 20 m/px |
| CRS | `EPSG:4326` |
| maxPixels Paso 1 | `1e9` |
| maxPixels Paso 2-3 | `1e13` |
| Filtro nubes | `CLOUDY_PIXEL_PERCENTAGE < 80` |
| Ventana gap-filling | ±30 días |
| Umbral LSP | 0.3 (30% amplitud relativa para SOS/EOS) |

---

## 8. Errores conocidos y soluciones aplicadas

| Error | Causa | Solución |
|---|---|---|
| `masked.divide is not a function` | `copyProperties()` devuelve `Element`, no `ee.Image` | Usar `.set('system:time_start', ...)` en lugar de `copyProperties` |
| Navegador GEE congelado | `.getInfo()` síncrono en loops (300+ llamadas) | Usar `.evaluate(callback)` asíncrono |
| `Cannot overwrite asset (Error code: 3)` | Asset ya existe de ejecución previa | Borrar asset previo o cambiar nombre con sufijo `_v2` |
| GPR calculado sobre tile completo | `.clip()` aplicado después del GPR | Aplicar `.clip(parcel)` antes de `makeGPR()` |
| `Expected ImageCollection, found IndexedFolder` | Carpetas creadas como Folder en vez de ImageCollection | Usar `folderToImageCollection()` con `ee.data.listAssets()` |
| `An asset with the ID already exists` | Task LSP previa creó asset vacío | Borrar asset previo desde panel Assets y hacer Run en la task |

---

## 9. Dependencias Python

Ver `environment.yml` en la raíz del repositorio.

Librerías principales: `earthengine-api`, `geemap`, `rasterio`, `numpy`, `pandas`, `matplotlib`, `geopandas`.
