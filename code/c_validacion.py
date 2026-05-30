# ============================================================
# Script C_v4: Validación — estilo IJRS/MDPI
#  
# ============================================================
import os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.features import geometry_mask
from shapely.geometry import Polygon, MultiPolygon, mapping
from shapely.ops import transform as shp_transform
from scipy import stats
from datetime import datetime
import pyproj
import warnings
warnings.filterwarnings("ignore")

BASE_DIR   = os.path.expanduser("~/GEE_Downloads")
OUTPUT_DIR = os.path.expanduser("~/analisis_gulupa/validacion_v3")
os.makedirs(OUTPUT_DIR, exist_ok=True)

N_POINTS  = 100
SEED      = 42
VARIABLES = ["LAI", "FVC", "laiCab"]
LABELS    = {
    "LAI":    "LAI (m² m⁻²)",
    "FVC":    "FVC (fracción)",
    "laiCab": "laiCab (µg cm⁻²)",
}

CMAP_VAR = {
    "LAI":    "YlGn",
    "FVC":    "GnBu",
    "laiCab": "YlOrBr",
}
CMAP_LABEL = {
    "LAI":    "LAI (m² m⁻²)",
    "FVC":    "FVC (fracción)",
    "laiCab": "laiCab (µg cm⁻²)",
}

# ── Geometría finca (AOI = Area_Total − Lote_Casa) ───────────
AREA_TOTAL_COORDS = [
    (-72.02911590018898, 4.323774005392068),
    (-72.02320164960223, 4.320341858927772),
    (-71.99426061433154, 4.337725140154523),
    (-71.99372417252856, 4.338110271794471),
    (-71.99303752702075, 4.338666572705499),
    (-71.99149257462817, 4.340207096166426),
    (-71.99072009843188, 4.340870376132908),
    (-71.98967671912509, 4.341426675010239),
    (-72.00730581274610, 4.349063111509516),
    (-72.00862576723136, 4.353153345576473),
    (-72.01213946104087, 4.356828060080524),
    (-72.01526010903413, 4.352303621556992),
    (-72.01894546422060, 4.346895815199964),
    (-72.02098296949191, 4.341000332940370),
    (-72.02443983146695, 4.334085615058952),
    (-72.02447201797513, 4.333775368489792),
    (-72.02456857749966, 4.333529310775441),
    (-72.02911590018898, 4.323774005392068),
]
LOTE_CASA_COORDS = [
    (-72.02527036942277, 4.327000978587103),
    (-72.02123566610653, 4.329945886052444),
    (-72.02102108938534, 4.329988678906903),
    (-72.02003403646786, 4.329753318177450),
    (-72.01848908407528, 4.329603543129691),
    (-72.01829596502621, 4.329582146691876),
    (-72.01656862242062, 4.328672797525615),
    (-72.01685830099423, 4.326244294986618),
    (-72.01699704976303, 4.325322370602150),
    (-72.01762468667252, 4.324659076998467),
    (-72.01774806828720, 4.324616283842795),
    (-72.01786608548386, 4.324659076998467),
    (-72.01965780110581, 4.325525637880559),
    (-72.01980264039261, 4.325616573224313),
    (-72.02000648827774, 4.325718206830860),
    (-72.02018351407273, 4.325793094742748),
    (-72.02030689568741, 4.325916124867642),
    (-72.02041954846604, 4.326097995450447),
    (-72.02046782822830, 4.326124741120704),
    (-72.02074141354782, 4.326172883324792),
    (-72.02128858418686, 4.326199628992403),
    (-72.02208251805527, 4.325456099080927),
    (-72.02527036942277, 4.327000978587103),
]
_POLY_TOTAL = Polygon(AREA_TOTAL_COORDS)
_POLY_CASA  = Polygon(LOTE_CASA_COORDS)
FINCA_POLY  = _POLY_TOTAL.difference(_POLY_CASA)
_minx, _miny, _maxx, _maxy = FINCA_POLY.bounds
_mx = (_maxx - _minx) * 0.04
_my = (_maxy - _miny) * 0.04
FINCA_XLIM = (_minx - _mx, _maxx + _mx)
FINCA_YLIM = (_miny - _my, _maxy + _my)

COLORS = {
    "GPR": "#1B4F72",
    "GF":  "#C0392B",
    "reg": "#E67E22",
    "pts": "#2E86AB",
}

# ── Estilo global artículo ───────────────────────────────────
plt.rcParams.update({
    "font.family":        "serif",
    "font.serif":         ["Times New Roman", "DejaVu Serif"],
    "font.size":          11,
    "axes.titlesize":     12,
    "axes.labelsize":     11,
    "xtick.labelsize":    10,
    "ytick.labelsize":    10,
    "legend.fontsize":    10,
    "legend.framealpha":  0.9,
    "legend.edgecolor":   "0.6",
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "axes.linewidth":     0.6,
    "xtick.major.width":  0.5,
    "ytick.major.width":  0.5,
    "xtick.direction":    "in",
    "ytick.direction":    "in",
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "grid.linewidth":     0.5,
    "grid.alpha":         0.35,
    "lines.linewidth":    1.5,
})


# ─── UTILIDADES ─────────────────────────────────────────────

def parse_date(filename):
    base = os.path.splitext(os.path.basename(filename))[0]
    m = re.search(r"(\d{4})-(\d{2})-(\d{2})", base)
    if m: return datetime(int(m.group(1)), int(m.group(2)), int(m.group(3)))
    m = re.search(r"(\d{8})", base)
    if m:
        s = m.group(1)
        return datetime(int(s[:4]), int(s[4:6]), int(s[6:8]))
    return None

def clean_array(arr, nodata=None):
    r = arr.astype(np.float64)
    if nodata is not None:
        r[r == nodata] = np.nan
    r[~np.isfinite(r)] = np.nan
    return r

def read_clean(filepath, band=1):
    with rasterio.open(filepath) as src:
        raw = src.read(band).astype(np.float64)
        return clean_array(raw, src.nodata), src.transform, src.crs, src.profile

def reproject_to_match(src_path, ref_path):
    with rasterio.open(ref_path) as ref:
        dst_crs, dst_tf = ref.crs, ref.transform
        dst_h, dst_w    = ref.height, ref.width
    with rasterio.open(src_path) as src:
        if src.crs == dst_crs and src.width == dst_w and src.height == dst_h:
            return clean_array(src.read(1).astype(np.float64), src.nodata)
        src_data = src.read(1).astype(np.float64)
        dst = np.full((dst_h, dst_w), np.nan, dtype=np.float64)
        reproject(
            source=src_data, destination=dst,
            src_transform=src.transform, src_crs=src.crs,
            dst_transform=dst_tf, dst_crs=dst_crs,
            resampling=Resampling.bilinear,
            src_nodata=src.nodata if src.nodata else -9999,
            dst_nodata=np.nan,
        )
    return clean_array(dst)

def best_reference_tif(folder, max_check=15):
    tifs = sorted([f for f in os.listdir(folder) if f.endswith(".tif")])
    sizes = {}
    for t in tifs[:max_check]:
        with rasterio.open(os.path.join(folder, t)) as src:
            key = (src.width, src.height)
            sizes[key] = sizes.get(key, 0) + 1
    modal_size = max(sizes, key=sizes.get)
    best_n, best_fp = 0, None
    for t in tifs[:max_check]:
        fp = os.path.join(folder, t)
        with rasterio.open(fp) as src:
            if (src.width, src.height) != modal_size:
                continue
            arr = clean_array(src.read(1).astype(np.float64), src.nodata)
        n_valid = int(np.sum(np.isfinite(arr) & (arr > 0)))
        if n_valid > best_n:
            best_n, best_fp = n_valid, fp
    print(f"    Ref. escogida : {os.path.basename(best_fp)}")
    print(f"    Píxeles válidos en ref: {best_n}")
    if best_n < N_POINTS:
        print(f"    ⚠  Solo {best_n} px válidos < {N_POINTS}.")
        print(f"       Aumenta max_check o revisa el TIF de referencia.")
    return best_fp


# ─── PASO 1: Generar puntos ─────────────────────────────────

def generate_points(ref_tif, n=100, seed=42):
    arr, transform, crs, _ = read_clean(ref_tif)
    valid = np.isfinite(arr) & (arr > 0)
    rows, cols = np.where(valid)
    rng = np.random.default_rng(seed)
    if len(rows) < n:
        print(f"  ⚠  Solo {len(rows)} px válidos → n={len(rows)}")
        n = len(rows)
    idx    = rng.choice(len(rows), size=n, replace=False)
    pts    = np.column_stack([rows[idx], cols[idx]])
    xs, ys = [], []
    for r, c in pts:
        x, y = rasterio.transform.xy(transform, int(r), int(c))
        xs.append(x); ys.append(y)
    return pts, np.array(xs), np.array(ys), crs


# ─── PASO 2: Extraer valores POR FECHA ──────────────────────

def extract_values(folder, points, ref_tif=None):
    tifs = sorted([f for f in os.listdir(folder) if f.endswith(".tif")])
    dates, matrix = [], []
    for tif in tifs:
        date = parse_date(tif)
        if date is None: continue
        fpath = os.path.join(folder, tif)
        data  = reproject_to_match(fpath, ref_tif) if ref_tif else \
                read_clean(fpath)[0]
        h, w  = data.shape
        vals  = [data[int(np.clip(r,0,h-1)), int(np.clip(c,0,w-1))]
                 for r, c in points]
        dates.append(date)
        matrix.append(vals)
    return dates, np.array(matrix, dtype=np.float64)


# ─── PASO 3: Métricas ───────────────────────────────────────

def metrics_per_date(dates_gpr, vals_gpr, dates_gf, vals_gf):
    records = []
    for i, d_gpr in enumerate(dates_gpr):
        best_j, best_diff = None, 999
        for j, d_gf in enumerate(dates_gf):
            diff = abs((d_gf - d_gpr).days)
            if diff < best_diff:
                best_diff, best_j = diff, j
        if best_diff > 5 or best_j is None:
            continue
        yt = vals_gpr[i];  yp = vals_gf[best_j]
        mask = np.isfinite(yt) & np.isfinite(yp) & (yt > 0) & (yp > 0)
        yt_m, yp_m = yt[mask], yp[mask]
        n_valid = len(yt_m)
        if n_valid < 5:
            continue
        residuals = yp_m - yt_m
        rmse   = float(np.sqrt(np.mean(residuals**2)))
        mae    = float(np.mean(np.abs(residuals)))
        bias   = float(np.mean(residuals))
        ss_res = float(np.sum(residuals**2))
        ss_tot = float(np.sum((yt_m - yt_m.mean())**2))
        r2     = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
        sl, ic, r_val, _, _ = stats.linregress(yt_m, yp_m) \
                               if n_valid >= 3 else (np.nan,)*5
        records.append({
            "date":      d_gpr,
            "n_valid":   n_valid,
            "RMSE":      round(rmse, 4),
            "MAE":       round(mae,  4),
            "R2":        round(float(r2), 4),
            "bias":      round(bias, 4),
            "slope":     round(float(sl), 4),
            "intercept": round(float(ic), 4),
            "pearson_r": round(float(r_val), 4),
        })
    return pd.DataFrame(records)

def metrics_global(df_per_date):
    cols = ["RMSE", "MAE", "R2", "bias", "pearson_r"]
    return df_per_date[cols].agg(["mean","std","median"]).round(4)


# ─── HELPERS MAPA ───────────────────────────────────────────

def _shapely_to_mpl_patch(geom, **kwargs):
    def _ring_codes(ring):
        v = list(ring.coords)
        c = [Path.MOVETO] + [Path.LINETO]*(len(v)-2) + [Path.CLOSEPOLY]
        return v, c
    all_v, all_c = [], []
    polys = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
    for poly in polys:
        v, c = _ring_codes(poly.exterior)
        all_v += v;  all_c += c
        for hole in poly.interiors:
            v, c = _ring_codes(hole)
            all_v += v;  all_c += c
    return PathPatch(Path(all_v, all_c), **kwargs)

def _mask_raster_to_finca(img, transform, crs_wkt):
    try:
        src_crs = pyproj.CRS.from_wkt(crs_wkt)
        wgs84   = pyproj.CRS.from_epsg(4326)
        if not src_crs.equals(wgs84):
            project    = pyproj.Transformer.from_crs(
                wgs84, src_crs, always_xy=True).transform
            finca_proj = shp_transform(project, FINCA_POLY)
        else:
            finca_proj = FINCA_POLY
        msk = geometry_mask(
            [mapping(finca_proj)],
            out_shape=img.shape,
            transform=transform,
            invert=True,
        )
        out = img.copy()
        out[~msk] = np.nan
        return out
    except Exception as e:
        print(f"    ⚠ _mask_raster_to_finca falló ({e})")
        return img


# ─── FIGURA A: Panel Scatter + Mapa ─────────────────────────

def _annotate_metrics(ax, m, n, loc="upper left"):
    txt = (
        f"$n$ = {n}\n"
        f"$R^2$ = {m['R2']:.3f}\n"
        f"RECM  = {m['RMSE']:.3f}\n"
        f"EAM   = {m['MAE']:.3f}\n"
        f"Sesgo = {m['bias']:.3f}\n"
        f"$r$   = {m['pearson_r']:.3f}\n"
        f"Pend. = {m['slope']:.3f}\n"
        f"Inter.= {m['intercept']:.3f}"
    )
    props = dict(boxstyle="round,pad=0.45", facecolor="white",
                 edgecolor="0.6", alpha=0.92, linewidth=0.6)
    x_pos = 0.03 if loc == "upper left" else 0.97
    ha    = "left" if loc == "upper left" else "right"
    ax.text(x_pos, 0.97, txt, transform=ax.transAxes,
            fontsize=9.5, verticalalignment="top",
            horizontalalignment=ha, bbox=props,
            family="serif")

def _draw_map_points(ax, xs, ys, ref_tif, var):
    with rasterio.open(ref_tif) as src:
        img   = src.read(1).astype(np.float64)
        img   = clean_array(img, src.nodata)
        tf    = src.transform
        crs_w = src.crs.to_wkt()
        left  = tf.c;  top   = tf.f
        right = left + tf.a * src.width
        bot   = top  + tf.e * src.height

    src_crs_obj  = pyproj.CRS.from_wkt(crs_w)
    is_geographic = src_crs_obj.is_geographic

    img_masked = _mask_raster_to_finca(img, tf, crs_w)
    p2,  p98   = np.nanpercentile(img_masked, 2), np.nanpercentile(img_masked, 98)
    img_clip   = np.clip(img_masked, p2, p98)

    cmap = CMAP_VAR.get(var, "YlGn")
    norm = Normalize(vmin=p2, vmax=p98)
    im   = ax.imshow(img_clip, extent=[left, right, bot, top],
                     origin="upper", cmap=cmap, norm=norm,
                     alpha=0.88, aspect="auto", zorder=1)

    divider = make_axes_locatable(ax)
    cax     = divider.append_axes("right", size="4%", pad=0.06)
    cb      = plt.colorbar(im, cax=cax)
    cb.set_label(CMAP_LABEL.get(var, var), fontsize=10, labelpad=4)
    cb.ax.tick_params(labelsize=9)
    cb.outline.set_linewidth(0.5)

    if is_geographic:
        ax.set_xlim(FINCA_XLIM);  ax.set_ylim(FINCA_YLIM)
        xs_plot, ys_plot = xs, ys
        finca_geom       = FINCA_POLY
        xlabel = "Longitud (°)";  ylabel = "Latitud (°)"
        fmt_x  = lambda v, _: f"{v:.4f}°"
        fmt_y  = lambda v, _: f"{v:.4f}°"
    else:
        wgs84   = pyproj.CRS.from_epsg(4326)
        project = pyproj.Transformer.from_crs(
            wgs84, src_crs_obj, always_xy=True).transform
        finca_geom = shp_transform(project, FINCA_POLY)
        fx0, fy0, fx1, fy1 = finca_geom.bounds
        mx = (fx1-fx0)*0.04;  my = (fy1-fy0)*0.04
        ax.set_xlim(fx0-mx, fx1+mx);  ax.set_ylim(fy0-my, fy1+my)
        xs_plot, ys_plot = xs, ys
        xlabel = "Este (m)";  ylabel = "Norte (m)"
        fmt_x  = lambda v, _: f"{v/1e3:.1f}k"
        fmt_y  = lambda v, _: f"{v/1e3:.1f}k"

    patch = _shapely_to_mpl_patch(
        finca_geom, facecolor="none",
        edgecolor="#333333", linewidth=0.9, zorder=4)
    ax.add_patch(patch)

    ax.scatter(xs_plot, ys_plot, s=18, c=COLORS["pts"],
               edgecolors="white", linewidths=0.4, zorder=5,
               label=f"Puntos de validación ($n$={len(xs_plot)})")

    ax.set_xlabel(xlabel, fontsize=10.5, labelpad=4)
    ax.set_ylabel(ylabel, fontsize=10.5, labelpad=4)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(fmt_x))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(fmt_y))
    ax.tick_params(labelsize=9.5)
    ax.set_facecolor("white")

    leg = ax.legend(fontsize=10, loc="lower right",
                    framealpha=0.9, edgecolor="0.6",
                    handletextpad=0.3, borderpad=0.4)
    leg.get_frame().set_linewidth(0.5)
    ax.set_title(f"Distribución espacial — {var}",
                 fontsize=11, pad=4)


def plot_panel_scatter_maps(all_data, ref_tifs):
    nrows = len(VARIABLES)
    fig   = plt.figure(figsize=(15, 4.6 * nrows))

    fig.suptitle(
        "Validación cruzada GPR vs. gap-filling (GF): diagramas de dispersión "
        "y distribución espacial de puntos de muestreo para LAI, FVC y laiCab",
        fontsize=12, fontweight="bold", y=0.995,
    )
    fig.subplots_adjust(
        left=0.08,  right=0.95,
        top=0.955,  bottom=0.05,
        hspace=0.30,
        wspace=0.35,
    )
    outer = gridspec.GridSpec(
        nrows, 2, figure=fig,
        width_ratios=[1, 1.15],
        hspace=0.30, wspace=0.35,
    )
    panel_letters = list("abcdef")
    letter_idx    = 0

    for row_i, var in enumerate(VARIABLES):
        if var not in all_data:
            continue

        d          = all_data[var]
        dates_gpr  = d["dates_gpr"];  vals_gpr  = d["vals_gpr"]
        dates_gf   = d["dates_gf"];   vals_gf   = d["vals_gf"]
        df_metrics = d["df_metrics"]; xs = d["xs"]; ys = d["ys"]

        ax_sc = fig.add_subplot(outer[row_i, 0])

        if df_metrics.empty:
            ax_sc.text(0.5, 0.5, "Sin datos",
                       transform=ax_sc.transAxes,
                       ha="center", va="center")
        else:
            best_row  = df_metrics.loc[df_metrics["R2"].idxmax()]
            best_date = best_row["date"]

            i_gpr = next((i for i, d_ in enumerate(dates_gpr)
                          if d_ == best_date), None)
            best_j, best_diff = None, 999
            for j, d_gf in enumerate(dates_gf):
                diff = abs((d_gf - best_date).days)
                if diff < best_diff:
                    best_diff, best_j = diff, j

            yt = vals_gpr[i_gpr];  yp = vals_gf[best_j]
            mask = np.isfinite(yt) & np.isfinite(yp) & (yt>0) & (yp>0)
            yt_m, yp_m = yt[mask], yp[mask]

            lims = [min(yt_m.min(), yp_m.min()) * 0.92,
                    max(yt_m.max(), yp_m.max()) * 1.08]

            ax_sc.scatter(yt_m, yp_m, s=32, alpha=0.60,
                          color=COLORS["pts"], edgecolors="white",
                          linewidths=0.35, zorder=3)
            ax_sc.plot(lims, lims, color="black", lw=1.0,
                       linestyle="--", zorder=2, label="Línea 1:1")
            if np.isfinite(best_row["slope"]):
                x_f = np.linspace(lims[0], lims[1], 300)
                ax_sc.plot(x_f,
                           best_row["slope"]*x_f + best_row["intercept"],
                           color=COLORS["reg"], lw=1.6, zorder=4,
                           label="Regresión")

            ax_sc.set_xlim(lims);  ax_sc.set_ylim(lims)
            ax_sc.set_xlabel(f"{var} — Referencia GPR (S2 directo)",
                             fontsize=10.5, labelpad=4)
            ax_sc.set_ylabel(f"{var} — Estimación GF (gap-filled)",
                             fontsize=10.5, labelpad=4)
            ax_sc.set_title(f"Validación dispersión — {var}",
                            fontsize=11, fontweight="bold", pad=4)
            ax_sc.legend(fontsize=10, loc="lower right",
                         framealpha=0.9, edgecolor="0.6",
                         handletextpad=0.4)
            ax_sc.grid(True, linestyle="--", alpha=0.35)
            _annotate_metrics(ax_sc, best_row, len(yt_m))

        ax_sc.text(-0.14, 1.04,
                   f"({panel_letters[letter_idx]})",
                   transform=ax_sc.transAxes,
                   fontsize=12, fontweight="bold", va="top")
        letter_idx += 1

        ax_mp = fig.add_subplot(outer[row_i, 1])
        ref   = ref_tifs.get(var)
        if ref and xs is not None and len(xs) > 0:
            _draw_map_points(ax_mp, xs, ys, ref, var)
        else:
            ax_mp.text(0.5, 0.5, "Mapa no disponible",
                       transform=ax_mp.transAxes,
                       ha="center", va="center")

        ax_mp.text(-0.12, 1.04,
                   f"({panel_letters[letter_idx]})",
                   transform=ax_mp.transAxes,
                   fontsize=12, fontweight="bold", va="top")
        letter_idx += 1

    out = os.path.join(OUTPUT_DIR, "FigA_scatter_maps.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Figura A guardada: {out}")


# ─── FIGURA 2: Serie temporal con percentiles ───────────────

def plot_temporal_percentiles(dates_gpr, vals_gpr,
                               dates_gf,  vals_gf, var):
    def pstats(matrix):
        m = matrix.copy()
        m[~np.isfinite(m) | (m <= 0)] = np.nan
        return (np.nanmedian(m, axis=1),
                np.nanpercentile(m, 25, axis=1),
                np.nanpercentile(m, 75, axis=1))

    med_gpr, p25_gpr, p75_gpr = pstats(vals_gpr)
    med_gf,  p25_gf,  p75_gf  = pstats(vals_gf)

    fig, ax = plt.subplots(figsize=(13, 4.8))
    fig.subplots_adjust(left=0.09, right=0.97, top=0.90, bottom=0.16)

    ax.fill_between(dates_gf, p25_gf, p75_gf,
                    alpha=0.20, color=COLORS["GF"],
                    label="GF IQR (P25–P75)")
    ax.plot(dates_gf, med_gf, color=COLORS["GF"], lw=2.0,
            label="GF mediana (gap-filled)")
    ax.fill_between(dates_gpr, p25_gpr, p75_gpr,
                    alpha=0.15, color=COLORS["GPR"])
    ax.plot(dates_gpr, med_gpr, "o-", color=COLORS["GPR"],
            markersize=5, lw=1.4,
            label="GPR mediana (S2 directo)")

    ax.set_xlabel("Fecha (2024)", fontsize=11, labelpad=6)
    ax.set_ylabel(LABELS[var], fontsize=11, labelpad=6)
    ax.set_title(
        f"Serie temporal — {var}  "
        f"(mediana ± IQR,  $n$ = {N_POINTS} puntos aleatorios)",
        fontsize=12, fontweight="bold", pad=8,
    )
    ax.legend(fontsize=10, loc="upper left",
              bbox_to_anchor=(0.01, 0.99),
              framealpha=0.9, edgecolor="0.6", ncol=1)
    ax.grid(True, linestyle="--", alpha=0.35)
    plt.xticks(rotation=30, ha="right", fontsize=10)

    out = os.path.join(OUTPUT_DIR, f"temporal_percentiles_{var}.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Serie temporal guardada: {out}")


# ─── FIGURA 3: Métricas por fecha ───────────────────────────

def plot_metrics_over_time(df_metrics, var):
    if df_metrics.empty: return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 6.5), sharex=True)
    fig.subplots_adjust(
        left=0.09, right=0.97,
        top=0.92,  bottom=0.12,
        hspace=0.18,
    )

    ax1.plot(df_metrics["date"], df_metrics["RMSE"],
             "o-", color=COLORS["GF"],
             markersize=4, lw=1.5, label="RECM")
    ax1.plot(df_metrics["date"], df_metrics["MAE"],
             "s--", color="#57A773",
             markersize=4, lw=1.2, label="EAM")
    ax1.set_ylabel(f"Error  ({LABELS[var]})", fontsize=11, labelpad=6)
    ax1.legend(fontsize=10, loc="upper right",
               framealpha=0.9, edgecolor="0.6")
    ax1.grid(True, linestyle="--", alpha=0.35)
    ax1.set_title(
        f"Métricas de validación por fecha — {var}  (n = 100 puntos / fecha)",
        fontsize=12, fontweight="bold", pad=8,
    )

    ax2.plot(df_metrics["date"], df_metrics["R2"],
             "o-", color=COLORS["GPR"],
             markersize=4, lw=1.5, label="$R^2$")
    ax2.axhline(0,   color="red",   linestyle="--", lw=1, alpha=0.6)
    ax2.axhline(0.7, color="green", linestyle="--", lw=1,
                alpha=0.6, label="$R^2 = 0{,}70$")
    ax2.set_ylabel("$R^2$", fontsize=11, labelpad=6)
    ax2.set_xlabel("Fecha (2024)", fontsize=11, labelpad=6)
    ax2.legend(fontsize=10, loc="lower right",
               framealpha=0.9, edgecolor="0.6")
    ax2.grid(True, linestyle="--", alpha=0.35)
    plt.xticks(rotation=30, ha="right", fontsize=10)

    out = os.path.join(OUTPUT_DIR, f"metricas_temporal_{var}.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Métricas temporales guardadas: {out}")


# ─── EJECUCIÓN PRINCIPAL ────────────────────────────────────

if __name__ == "__main__":
    print("=" * 60)
    print("  Validación v4 — n=100 | IJRS/MDPI | Times New Roman")
    print("=" * 60)

    summary_global = []
    all_data = {}
    ref_tifs = {}

    for var in VARIABLES:
        print(f"\n{'─'*55}\n  Variable: {var}\n{'─'*55}")

        folder_gpr = os.path.join(BASE_DIR, f"GPR_{var}_2024")
        folder_gf  = os.path.join(BASE_DIR, f"GF_{var}_2024")

        if not os.path.isdir(folder_gpr) or not os.path.isdir(folder_gf):
            print(f"  ⚠  Carpetas no encontradas"); continue

        ref_tif = best_reference_tif(folder_gpr)
        ref_tifs[var] = ref_tif

        points, xs, ys, crs = generate_points(ref_tif, n=N_POINTS, seed=SEED)
        pd.DataFrame({"punto_id": range(len(points)),
                      "fila": points[:,0], "col": points[:,1],
                      "lon": xs, "lat": ys})\
          .to_csv(os.path.join(OUTPUT_DIR, f"puntos_100_{var}.csv"),
                  index=False)

        print("  Extrayendo GPR...")
        dates_gpr, vals_gpr = extract_values(folder_gpr, points)
        print("  Extrayendo GF (reproyectando)...")
        dates_gf, vals_gf = extract_values(folder_gf, points,
                                            ref_tif=ref_tif)

        df_m = metrics_per_date(dates_gpr, vals_gpr, dates_gf, vals_gf)
        df_m.to_csv(os.path.join(OUTPUT_DIR,
                                  f"metricas_por_fecha_{var}.csv"),
                    index=False)

        gsum = metrics_global(df_m)
        print(f"\n  Resumen global {var}:")
        print(gsum.to_string())
        summary_global.append({
            "variable": var,
            **{f"{c}_media": gsum.loc["mean", c] for c in gsum.columns},
            **{f"{c}_de":    gsum.loc["std",  c] for c in gsum.columns},
        })

        all_data[var] = dict(
            dates_gpr=dates_gpr, vals_gpr=vals_gpr,
            dates_gf=dates_gf,   vals_gf=vals_gf,
            df_metrics=df_m,     xs=xs, ys=ys,
        )

        plot_temporal_percentiles(dates_gpr, vals_gpr,
                                   dates_gf,  vals_gf, var)
        plot_metrics_over_time(df_m, var)

    if all_data:
        plot_panel_scatter_maps(all_data, ref_tifs)

    pd.DataFrame(summary_global)\
      .to_csv(os.path.join(OUTPUT_DIR, "resumen_global_v4.csv"),
              index=False)
    print(f"\n✅ Todo guardado en: {OUTPUT_DIR}")
    
