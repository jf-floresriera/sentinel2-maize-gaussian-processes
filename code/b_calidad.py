# ============================================================
# Script b: Caracterización de calidad de imágenes descargadas
# ============================================================
import os, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rasterio
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")

BASE_DIR   = os.path.expanduser("~/GEE_Downloads")
OUTPUT_DIR = os.path.expanduser("~/analisis_gulupa/calidad")
os.makedirs(OUTPUT_DIR, exist_ok=True)

VARIABLES = ["LAI", "FVC", "laiCab"]
STEPS     = ["GPR", "GF"]

def parse_date(filename):
    base = os.path.splitext(os.path.basename(filename))[0]
    m = re.search(r"(\d{4})-(\d{2})-(\d{2})", base)
    if m: return datetime(int(m.group(1)), int(m.group(2)), int(m.group(3)))
    m = re.search(r"(\d{8})", base)
    if m:
        s = m.group(1)
        return datetime(int(s[:4]), int(s[4:6]), int(s[6:8]))
    return None
def quality_metrics(filepath):
    """Calcula métricas de calidad para un GeoTIFF."""
    with rasterio.open(filepath) as src:
        data   = src.read(1).astype(np.float32)
        nodata = src.nodata
        meta   = {
            "width":      src.width,
            "height":     src.height,
            "crs":        str(src.crs),
            "resolution": abs(src.transform[0]),
            "n_bands":    src.count,
        }
    if nodata is not None:
        data[data == nodata] = np.nan
    data[np.isneginf(data)] = np.nan   # ← LÍNEA NUEVA: elimina -inf de GEE
    data[np.isposinf(data)] = np.nan   # ← LÍNEA NUEVA: elimina +inf de GEE
    
    total  = data.size
    valids = int(np.sum(~np.isnan(data)))
    nans   = total - valids

    return {
        **meta,
        "total_pixels":   total,
        "valid_pixels":   valids,
        "nodata_pixels":  nans,
        "coverage_%":     round(valids / total * 100, 2),  # calidad clave
        "mean":           round(float(np.nanmean(data)), 4),
        "std":            round(float(np.nanstd(data)), 4),
        "min":            round(float(np.nanmin(data)), 4) if valids > 0 else np.nan,
        "max":            round(float(np.nanmax(data)), 4) if valids > 0 else np.nan,
        "p25":            round(float(np.nanpercentile(data, 25)), 4) if valids > 0 else np.nan,
        "median":         round(float(np.nanpercentile(data, 50)), 4) if valids > 0 else np.nan,
        "p75":            round(float(np.nanpercentile(data, 75)), 4) if valids > 0 else np.nan,
        "cv_%":           round(float(np.nanstd(data) / np.nanmean(data) * 100), 2)
                          if valids > 0 and np.nanmean(data) != 0 else np.nan,
    }

# ─── ANÁLISIS DE CALIDAD POR COLECCIÓN ───────────────────────
all_records = []

for var in VARIABLES:
    for step in STEPS:
        folder = os.path.join(BASE_DIR, f"{step}_{var}_2024")
        if not os.path.isdir(folder):
            continue
        tifs = sorted([f for f in os.listdir(folder) if f.endswith(".tif")])
        print(f"\n📊 {step}_{var}_2024  →  {len(tifs)} imágenes")

        coverages, dates = [], []
        for tif in tifs:
            fpath = os.path.join(folder, tif)
            date  = parse_date(tif)
            try:
                metrics = quality_metrics(fpath)
                record  = {"variable": var, "step": step,
                           "date": date, "filename": tif, **metrics}
                all_records.append(record)
                coverages.append(metrics["coverage_%"])
                if date: dates.append(date)
                print(f"  {tif:<35} coverage={metrics['coverage_%']:>6.1f}%  "
                      f"mean={metrics['mean']:>8.4f}  std={metrics['std']:>7.4f}")
            except Exception as e:
                print(f"  ❌ Error: {tif} → {e}")

        # ── Gráfico de cobertura temporal por colección ──────
        if dates and coverages:
            fig, ax = plt.subplots(figsize=(11, 3.5))
            ax.bar(dates, coverages, width=5, color="#2E86AB", alpha=0.75,
                   edgecolor="white")
            ax.axhline(80, color="red",    linestyle="--", linewidth=1,
                       label="Umbral 80%")
            ax.axhline(95, color="green",  linestyle="--", linewidth=1,
                       label="Umbral 95%")
            ax.set_ylim(0, 105)
            ax.set_ylabel("Cobertura válida (%)")
            ax.set_xlabel("Fecha")
            ax.set_title(f"Calidad de imágenes — {step}_{var}_2024",
                         fontweight="bold")
            ax.legend(fontsize=9)
            ax.grid(True, axis="y", linestyle="--", alpha=0.4)
            plt.xticks(rotation=30, ha="right")
            plt.tight_layout()
            out = os.path.join(OUTPUT_DIR,
                               f"calidad_cobertura_{step}_{var}.png")
            plt.savefig(out, dpi=180, bbox_inches="tight")
            plt.close()
            print(f"  💾 {out}")

# ─── LSP: Calidad de imágenes multi-banda ────────────────────
LSP_BANDS_ORDER = ["sos", "eos", "pos", "los"]
for var in VARIABLES:
    lsp_file = os.path.join(BASE_DIR, f"LSP_{var}_2024.tif")
    if not os.path.exists(lsp_file):
        continue
    print(f"\n🗺  LSP_{var}_2024.tif")
    with rasterio.open(lsp_file) as src:
        for i in range(1, src.count + 1):
            band_name = src.descriptions[i-1] or LSP_BANDS_ORDER[i-1] \
                        if i <= len(LSP_BANDS_ORDER) else f"band_{i}"
            arr    = src.read(i).astype(np.float32)
            nodata = src.nodata
            if nodata is not None:
                arr[arr == nodata] = np.nan
            cov = np.sum(~np.isnan(arr)) / arr.size * 100
            print(f"  Banda {band_name:<6} → coverage={cov:.1f}%  "
                  f"mean={np.nanmean(arr):.2f}  "
                  f"range=[{np.nanmin(arr):.1f}, {np.nanmax(arr):.1f}]")

# ─── CSV RESUMEN COMPLETO ─────────────────────────────────────
if all_records:
    df = pd.DataFrame(all_records)
    df["date"] = pd.to_datetime(df["date"])
    df = df.sort_values(["variable", "step", "date"])
    csv_out = os.path.join(OUTPUT_DIR, "resumen_calidad_imagenes.csv")
    df.to_csv(csv_out, index=False)
    print(f"\n✅ CSV guardado: {csv_out}")

    # ── Figura resumen: boxplot cobertura por colección ──────
    fig, ax = plt.subplots(figsize=(12, 5))
    groups  = [f"{step}_{var}" for var in VARIABLES for step in STEPS]
    data_bp = []
    labels_bp = []
    for var in VARIABLES:
        for step in STEPS:
            sub = df[(df["variable"]==var) & (df["step"]==step)]["coverage_%"]
            if len(sub):
                data_bp.append(sub.values)
                labels_bp.append(f"{step}\n{var}")

    colors_bp = ["#2E86AB","#57A773"] * len(VARIABLES)
    bp = ax.boxplot(data_bp, patch_artist=True,
                    medianprops=dict(color="black", linewidth=2))
    for patch, col in zip(bp["boxes"], colors_bp):
        patch.set_facecolor(col)
        patch.set_alpha(0.7)
    ax.set_xticklabels(labels_bp, fontsize=9)
    ax.axhline(80, color="red",  linestyle="--", linewidth=1, label="80%")
    ax.axhline(95, color="green",linestyle="--", linewidth=1, label="95%")
    ax.set_ylabel("Cobertura válida (%)")
    ax.set_title("Calidad general — GPR vs GF por variable", fontweight="bold")
    ax.legend(); ax.grid(True, axis="y", alpha=0.4)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "calidad_resumen_boxplot.png"),
                dpi=180, bbox_inches="tight")
    plt.close()
    print("✅ Boxplot resumen guardado")
    