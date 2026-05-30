# ============================================================
# Figuras faltantes para artículo — estilo IJRS/MDPI
#   1. Mapa de ubicación (S. América → Colombia → finca)
#   2. Mapas GPR vs GF en etapas fenológicas clave (Ciclo B)
#   3. Figura comparativa estilo Figure 6 (Salinero-Delgado 2022)
#   4. Integración fenología + imágenes satelitales
# ============================================================
import os, re, json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe
from matplotlib.colors import Normalize
from matplotlib.path import Path
from matplotlib.patches import PathPatch, FancyArrowPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.font_manager as fm
import rasterio
from rasterio.warp import reproject, Resampling
import rasterio.features
from shapely.geometry import Polygon, mapping
import geopandas as gpd
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")

# ── Intentar importar cartopy ────────────────────────────────
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("⚠  cartopy no instalado — mapa de ubicación usará matplotlib simple")

# ── Intentar importar contextily ─────────────────────────────
try:
    import contextily as ctx
    HAS_CTX = True
except ImportError:
    HAS_CTX = False
    print("⚠  contextily no instalado — fondo satelital no disponible")


# ─── ESTILO GLOBAL ARTÍCULO ──────────────────────────────────
plt.rcParams.update({
    "font.family":        "serif",
    "font.serif":         ["Times New Roman", "DejaVu Serif"],
    "font.size":          9,
    "axes.titlesize":     10,
    "axes.labelsize":     9,
    "xtick.labelsize":    8,
    "ytick.labelsize":    8,
    "legend.fontsize":    8,
    "legend.framealpha":  0.9,
    "legend.edgecolor":   "0.6",
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "axes.linewidth":     0.6,
    "xtick.major.width":  0.5,
    "ytick.major.width":  0.5,
    "xtick.direction":    "in",
    "ytick.direction":    "in",
})

BASE_DIR   = os.path.expanduser("~/GEE_Downloads")
OUTPUT_DIR = os.path.expanduser("~/analisis_gulupa/figuras_articulo")
os.makedirs(OUTPUT_DIR, exist_ok=True)

VARIABLES = ["LAI", "FVC", "laiCab"]
CMAPS     = {"LAI": "YlGn", "FVC": "YlGn", "laiCab": "RdYlGn"}
VRANGES   = {"LAI": (0, 6), "FVC": (0, 1), "laiCab": (0, 2.5)}
UNITS     = {"LAI": "m² m⁻²", "FVC": "fracción", "laiCab": "g m⁻²"}

# ─── POLÍGONO REAL DE LA FINCA (WGS84) ───────────────────────
AREA_TOTAL_COORDS = [
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
    [-72.02911590018898, 4.323774005392068],
]
LOTE_CASA_COORDS = [
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
    [-72.02527036942277, 4.327000978587103],
]

_poly_total = Polygon(AREA_TOTAL_COORDS)
_poly_casa  = Polygon(LOTE_CASA_COORDS)
FINCA_POLY  = _poly_total.difference(_poly_casa)
FINCA_GDF   = gpd.GeoDataFrame(geometry=[FINCA_POLY], crs="EPSG:4326")
FINCA_LON, FINCA_LAT = FINCA_POLY.centroid.x, FINCA_POLY.centroid.y

# ─── ETAPAS FENOLÓGICAS ───────────────────────────────────────
PHENOSTAGES = {
    "VE":  ("Emergencia",           7,  113, 240),
    "V4":  ("4ª hoja",             20,  126, 253),
    "VT":  ("Panoja",              70,  176, 303),
    "R1":  ("Espigado",            75,  181, 308),
    "R3":  ("Grano lechoso",       95,  201, 328),
    "R6":  ("Madurez fisiológica", 120, 226, 353),
}

TARGET_DOYS_B = [253, 303, 328]
STAGE_NAMES_B = [
    "V4 — 4ª hoja\n(Veg.)",
    "VT — Panoja\n(Pico)",
    "R3 — Grano lechoso\n(Reprod.)",
]


# ─── UTILIDADES ──────────────────────────────────────────────
def parse_date(filename):
    base = os.path.splitext(os.path.basename(filename))[0]
    m = re.search(r"(\d{4})-(\d{2})-(\d{2})", base)
    if m: return datetime(int(m.group(1)), int(m.group(2)), int(m.group(3)))
    m = re.search(r"(\d{8})", base)
    if m:
        s = m.group(1)
        return datetime(int(s[:4]), int(s[4:6]), int(s[6:8]))
    return None


def clean(arr, nodata=None):
    r = arr.astype(np.float32)
    if nodata is not None: r[r == nodata] = np.nan
    r[~np.isfinite(r)] = np.nan
    return r


def modal_shape(folder, n=20):
    tifs = sorted([f for f in os.listdir(folder) if f.endswith(".tif")])[:n]
    s = {}
    for t in tifs:
        with rasterio.open(os.path.join(folder, t)) as src:
            k = (src.width, src.height)
            s[k] = s.get(k, 0) + 1
    return max(s, key=s.get)


def load_collection(folder):
    msh   = modal_shape(folder)
    tifs  = sorted([f for f in os.listdir(folder) if f.endswith(".tif")])
    dates, arrays, transforms = [], [], []
    for tif in tifs:
        date = parse_date(tif)
        if date is None: continue
        fp = os.path.join(folder, tif)
        with rasterio.open(fp) as src:
            if (src.width, src.height) != msh: continue
            arr = clean(src.read(1).astype(np.float32), src.nodata)
            tf  = src.transform
        dates.append(date)
        arrays.append(arr)
        transforms.append(tf)
    pairs = sorted(zip(dates, arrays, transforms), key=lambda x: x[0])
    if pairs:
        dates, arrays, transforms = zip(*pairs)
    return list(dates), list(arrays), list(transforms)


def reproject_to_ref(src_path, ref_path):
    with rasterio.open(ref_path) as ref:
        dst_crs, dst_tf = ref.crs, ref.transform
        dst_h, dst_w    = ref.height, ref.width
    with rasterio.open(src_path) as src:
        if src.crs == dst_crs and (src.width, src.height) == (dst_w, dst_h):
            return clean(src.read(1).astype(np.float32), src.nodata)
        sdata = src.read(1).astype(np.float32)
        dst   = np.full((dst_h, dst_w), np.nan, dtype=np.float32)
        reproject(source=sdata, destination=dst,
                  src_transform=src.transform, src_crs=src.crs,
                  dst_transform=dst_tf, dst_crs=dst_crs,
                  resampling=Resampling.bilinear,
                  src_nodata=src.nodata or -9999, dst_nodata=np.nan)
    return clean(dst)


def get_aoi(arr, pad=4):
    v    = np.isfinite(arr) & (arr > 0)
    rows = np.where(v.any(axis=1))[0]
    cols = np.where(v.any(axis=0))[0]
    if len(rows) == 0: return 0, arr.shape[0], 0, arr.shape[1]
    return (max(0, rows[0]-pad), min(arr.shape[0], rows[-1]+pad),
            max(0, cols[0]-pad), min(arr.shape[1], cols[-1]+pad))


def crop(arr, ext):
    return arr[ext[0]:ext[1], ext[2]:ext[3]]


def panel_label(ax, ltr, fs=9):
    ax.text(0.02, 0.97, f"({ltr})", transform=ax.transAxes,
            fontsize=fs, fontweight="bold", va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.15", fc="white",
                      ec="none", alpha=0.85))


def colorbar_h(fig, ax, im, label, nticks=5):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="6%", pad=0.10)
    cb  = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label(label, fontsize=7, labelpad=1)
    cb.ax.tick_params(labelsize=7)
    cb.locator = mticker.MaxNLocator(nbins=nticks)
    cb.update_ticks()
    return cb


def find_closest_date(target_doy, dates):
    doys = [(abs(d.timetuple().tm_yday - target_doy), i)
            for i, d in enumerate(dates)]
    return min(doys, key=lambda x: x[0])[1]


def mask_array_to_finca(arr, src_transform, src_crs):
    """Enmascara a NaN todos los píxeles fuera de FINCA_POLY."""
    h, w = arr.shape
    try:
        epsg = int(src_crs.to_epsg())
        finca_repr = FINCA_GDF.to_crs(epsg=epsg)
    except Exception:
        finca_repr = FINCA_GDF
    mask_arr = rasterio.features.geometry_mask(
        [mapping(geom) for geom in finca_repr.geometry],
        out_shape=(h, w),
        transform=src_transform,
        invert=True,
        all_touched=True,
    )
    result = arr.copy()
    result[~mask_arr] = np.nan
    return result


def add_north_arrow(ax, x=0.92, y=0.15, size=0.06, color="black"):
    """Flecha norte geográfico en coordenadas de ejes."""
    ax.annotate(
        "", xy=(x, y + size), xytext=(x, y),
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="-|>", color=color,
                        lw=1.5, mutation_scale=12),
    )
    ax.text(x, y + size + 0.03, "N", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=8,
            fontweight="bold", color=color)


def add_scale_bar_deg(ax, lon0, lat0, length_deg, label, color="black"):
    """Barra de escala simple en coordenadas geográficas."""
    ax.plot([lon0, lon0 + length_deg], [lat0, lat0],
            color=color, lw=2, transform=ax.transData, zorder=10)
    ax.text(lon0 + length_deg / 2, lat0 - 0.002, label,
            ha="center", va="top", fontsize=7, color=color,
            transform=ax.transData, zorder=10)


LETTERS = list("abcdefghijklmnopqrstuvwxyz")


# ════════════════════════════════════════════════════════════
# FIGURA A: Mapa de Ubicación — 3 paneles:
#   (a) América del Sur   (b) Colombia   (c) Finca satelital
# ════════════════════════════════════════════════════════════
def fig_location_map():
    fig = plt.figure(figsize=(12, 4.5))
    gs  = gridspec.GridSpec(1, 3, figure=fig,
                            wspace=0.08, left=0.02, right=0.98,
                            top=0.88, bottom=0.05)

    subtitles = ["(a) América del Sur",
                 "(b) Colombia",
                 "(c) Sitio de estudio (Meta, Colombia)"]

    # ── Panel (a): América del Sur ────────────────────────────
    if HAS_CARTOPY:
        prj = ccrs.PlateCarree()
        ax_sa = fig.add_subplot(gs[0], projection=prj)
        ax_sa.set_extent([-85, -32, -58, 15], crs=prj)
        ax_sa.add_feature(cfeature.LAND,      facecolor="#e8e0d0",
                          edgecolor="0.5",    linewidth=0.3)
        ax_sa.add_feature(cfeature.OCEAN,     facecolor="#c9dff0")
        ax_sa.add_feature(cfeature.BORDERS,   linewidth=0.4, edgecolor="0.4")
        ax_sa.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor="0.4")
        # Resaltar Colombia
        ax_sa.add_feature(
            cfeature.NaturalEarthFeature(
                "cultural", "admin_0_countries", "50m",
                facecolor="#f4b942", edgecolor="0.4", linewidth=0.5),
        )
        # Rectángulo Colombia
        col_ext = [-82, -66, -5, 13]
        rect_col = mpatches.Rectangle(
            (col_ext[0], col_ext[2]),
            col_ext[1]-col_ext[0], col_ext[3]-col_ext[2],
            linewidth=1.2, edgecolor="red", facecolor="none",
            transform=prj, zorder=5)
        ax_sa.add_patch(rect_col)
        # Punto finca
        ax_sa.plot(FINCA_LON, FINCA_LAT, "r*", ms=6,
                   transform=prj, zorder=6)
        gl = ax_sa.gridlines(draw_labels=True, linewidth=0.3,
                              color="gray", alpha=0.4,
                              x_inline=False, y_inline=False)
        gl.top_labels = False; gl.right_labels = False
        gl.xlabel_style = {"size": 6}; gl.ylabel_style = {"size": 6}
    else:
        ax_sa = fig.add_subplot(gs[0])
        ax_sa.set_facecolor("#c9dff0")
        ax_sa.set_xlim(-85, -32); ax_sa.set_ylim(-58, 15)
        ax_sa.text(0.5, 0.5, "América del Sur\n(instalar cartopy)",
                   transform=ax_sa.transAxes, ha="center")

    ax_sa.set_title(subtitles[0], fontsize=8.5, fontweight="bold", pad=3)

    # ── Panel (b): Colombia ───────────────────────────────────
    if HAS_CARTOPY:
        ax_col = fig.add_subplot(gs[1], projection=prj)
        ax_col.set_extent([-82, -66, -5, 13], crs=prj)
        ax_col.add_feature(cfeature.LAND,      facecolor="#e8e0d0",
                           edgecolor="0.5",    linewidth=0.3)
        ax_col.add_feature(cfeature.OCEAN,     facecolor="#c9dff0")
        ax_col.add_feature(cfeature.BORDERS,   linewidth=0.5, edgecolor="0.4")
        ax_col.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor="0.4")
        ax_col.add_feature(cfeature.RIVERS,    linewidth=0.4,
                           edgecolor="#7ecbdb", alpha=0.6)
        # Departamento Meta aproximado
        meta_lon, meta_lat = FINCA_LON, FINCA_LAT
        buf = 1.2
        meta_rect = mpatches.Rectangle(
            (meta_lon - buf, meta_lat - buf),
            buf * 2, buf * 2,
            linewidth=1.2, edgecolor="red", facecolor="#f4b942",
            alpha=0.35, transform=prj, zorder=4)
        ax_col.add_patch(meta_rect)
        ax_col.plot(FINCA_LON, FINCA_LAT, "r*", ms=7,
                    transform=prj, zorder=6)
        gl2 = ax_col.gridlines(draw_labels=True, linewidth=0.3,
                               color="gray", alpha=0.4,
                               x_inline=False, y_inline=False)
        gl2.top_labels = False; gl2.right_labels = False
        gl2.xlabel_style = {"size": 6}; gl2.ylabel_style = {"size": 6}
        # Etiqueta "Meta"
        ax_col.text(FINCA_LON + 0.3, FINCA_LAT - 0.8, "Meta",
                    transform=prj, fontsize=7.5, color="#8B0000",
                    fontweight="bold", zorder=7)
    else:
        ax_col = fig.add_subplot(gs[1])
        ax_col.set_facecolor("#c9dff0")
        ax_col.text(0.5, 0.5, "Colombia\n(instalar cartopy)",
                    transform=ax_col.transAxes, ha="center")

    ax_col.set_title(subtitles[1], fontsize=8.5, fontweight="bold", pad=3)

    # ── Panel (c): Finca con imagen satelital ─────────────────
    ax_d = fig.add_subplot(gs[2])

    finca_web  = FINCA_GDF.to_crs(epsg=3857)
    area_gdf   = gpd.GeoDataFrame(
        geometry=[_poly_total], crs="EPSG:4326").to_crs(epsg=3857)
    casa_gdf   = gpd.GeoDataFrame(
        geometry=[_poly_casa],  crs="EPSG:4326").to_crs(epsg=3857)

    # Establecer extent antes del basemap
    bounds_w   = finca_web.total_bounds        # (minx,miny,maxx,maxy)
    pad_m      = 350                           # metros de margen
    ax_d.set_xlim(bounds_w[0] - pad_m, bounds_w[2] + pad_m)
    ax_d.set_ylim(bounds_w[1] - pad_m, bounds_w[3] + pad_m)

    if HAS_CTX:
        try:
            ctx.add_basemap(
                ax_d,
                crs=finca_web.crs.to_string(),
                source=ctx.providers.Esri.WorldImagery,
                zoom=16,
                attribution=False,
            )
        except Exception:
            try:
                ctx.add_basemap(
                    ax_d,
                    crs=finca_web.crs.to_string(),
                    source=ctx.providers.OpenStreetMap.Mapnik,
                    zoom=15,
                    attribution=False,
                )
            except Exception:
                ax_d.set_facecolor("#4a7c59")

    # Linderos sobre el mosaico
    area_gdf.boundary.plot(ax=ax_d, color="red",    linewidth=1.8, zorder=6)
    finca_web.boundary.plot(ax=ax_d, color="yellow", linewidth=1.2,
                            linestyle="--", zorder=7)
    casa_gdf.plot(ax=ax_d,           color="white",  alpha=0.50,   zorder=5)
    casa_gdf.boundary.plot(ax=ax_d,  color="white",  linewidth=0.9, zorder=6)

    # ── Coordenadas en los ejes (Web Mercator → grados aprox.) ─
    from pyproj import Transformer
    tr = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)

    xlim = ax_d.get_xlim()
    ylim = ax_d.get_ylim()

    # 4 ticks en X y 4 en Y (Web Mercator) → convertir a grados
    xticks_m = np.linspace(xlim[0], xlim[1], 4)
    yticks_m = np.linspace(ylim[0], ylim[1], 4)
    xlabels  = [f"{tr.transform(x, ylim[0])[0]:.4f}°O" for x in xticks_m]
    ylabels  = [f"{tr.transform(xlim[0], y)[1]:.4f}°N" for y in yticks_m]

    ax_d.set_xticks(xticks_m)
    ax_d.set_xticklabels(xlabels, fontsize=6.0, rotation=30, ha="right",
                         color="black")
    ax_d.set_yticks(yticks_m)
    ax_d.set_yticklabels(ylabels, fontsize=6.0, color="black")
    ax_d.tick_params(axis="both", direction="out", length=3,
                     width=0.5, color="black", labelcolor="black")
    for spine in ax_d.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(0.7)

    # ── Norte geográfico ──────────────────────────────────────
    # Flecha dentro del panel (esquina superior derecha)
    ax_d.annotate(
        "", xy=(0.93, 0.93), xytext=(0.93, 0.80),
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="-|>", color="white",
                        lw=1.8, mutation_scale=14),
        zorder=15,
    )
    ax_d.text(0.93, 0.96, "N", transform=ax_d.transAxes,
              ha="center", va="bottom", fontsize=9,
              fontweight="bold", color="white", zorder=15,
              path_effects=[pe.withStroke(linewidth=1.5, foreground="black")])

    # ── Leyenda fuera del mapa (debajo) ──────────────────────
    legend_elements = [
        mpatches.Patch(facecolor="none", edgecolor="red",
                       linewidth=1.8, label="Perímetro total"),
        mpatches.Patch(facecolor="none", edgecolor="yellow",
                       linewidth=1.2, linestyle="--", label="Lote casa"),
        mpatches.Patch(facecolor="white", edgecolor="white",
                       label="Área maiz"),
    ]
    ax_d.legend(
        handles=legend_elements,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),   # debajo del panel
        ncol=3,
        fontsize=7,
        framealpha=0.95,
        edgecolor="0.5",
        facecolor="white",
        labelcolor="black",
    )

    ax_d.set_title(subtitles[2], fontsize=8.5, fontweight="bold", pad=3)

    fig.suptitle(
        "Ubicación del área de estudio — departamento del Meta, Colombia",
        fontsize=11, fontweight="bold"
    )
    out = os.path.join(OUTPUT_DIR, "FigLoc_location_map.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  ✅ {out}")


# ════════════════════════════════════════════════════════════
# FIGURA B: GPR → GF estilo Figure 6 — CICLO B
# ════════════════════════════════════════════════════════════
def fig_gpr_gf_comparison():
    for var in VARIABLES:
        folder_gpr = os.path.join(BASE_DIR, f"GPR_{var}_2024")
        folder_gf  = os.path.join(BASE_DIR, f"GF_{var}_2024")
        if not os.path.isdir(folder_gpr): continue

        dates_gpr, arr_gpr, _ = load_collection(folder_gpr)
        dates_gf,  arr_gf,  _ = load_collection(folder_gf)
        if not dates_gpr: continue

        gf_tifs = sorted([f for f in os.listdir(folder_gf)
                          if f.endswith(".tif")])
        msh = modal_shape(folder_gf)
        ref_gf = next(
            os.path.join(folder_gf, t) for t in gf_tifs
            if (lambda s: (s.width, s.height) == msh)(
                rasterio.open(os.path.join(folder_gf, t)))
        )

        extent = get_aoi(arr_gf[0], pad=5)
        vmin, vmax = VRANGES[var]
        cmap       = CMAPS[var]
        n_cols     = len(TARGET_DOYS_B)

        fig, axes = plt.subplots(
            3, n_cols,
            figsize=(3.5 * n_cols, 9.5),
            gridspec_kw={"hspace": 0.38, "wspace": 0.10}
        )

        row_titles = [
            f"{var} — Predicción GPR\n(observación directa S2)",
            f"{var} — Relleno de gaps GPR\n(reconstrucción temporal)",
            "Diferencia  GF − GPR",
        ]
        letter_idx = 0

        for col, (target_doy, stage_name) in enumerate(
                zip(TARGET_DOYS_B, STAGE_NAMES_B)):

            idx_gpr = find_closest_date(target_doy, dates_gpr)
            idx_gf  = find_closest_date(target_doy, dates_gf)
            date_gpr = dates_gpr[idx_gpr]
            date_gf  = dates_gf[idx_gf]

            gpr_fp = os.path.join(
                folder_gpr,
                sorted([f for f in os.listdir(folder_gpr)
                        if f.endswith(".tif")])[idx_gpr])
            gpr_arr = reproject_to_ref(gpr_fp, ref_gf)
            gf_arr  = arr_gf[idx_gf]

            gpr_c = crop(gpr_arr, extent)
            gf_c  = crop(gf_arr,  extent)
            dif_c = gf_c - gpr_c
            gpr_c = np.where(gpr_c > 0, gpr_c, np.nan)
            gf_c  = np.where(gf_c  > 0, gf_c,  np.nan)

            for row_idx, (data, title_txt) in enumerate([
                (gpr_c, f"{stage_name}\n{date_gpr.strftime('%d %b %Y')}"),
                (gf_c,  f"Relleno de gaps\n{date_gf.strftime('%d %b %Y')}"),
                (dif_c, "GF − GPR"),
            ]):
                ax = axes[row_idx, col]
                if row_idx < 2:
                    _cmap = plt.cm.get_cmap(cmap).copy()
                    _cmap.set_bad("white")
                    im = ax.imshow(np.ma.masked_invalid(data),
                                   cmap=_cmap, vmin=vmin, vmax=vmax,
                                   interpolation="nearest", aspect="equal")
                else:
                    d_abs = np.nanmax(np.abs(dif_c)) if np.any(np.isfinite(dif_c)) else 1
                    _cmap = plt.cm.get_cmap("RdBu_r").copy()
                    _cmap.set_bad("white")
                    im = ax.imshow(np.ma.masked_invalid(data),
                                   cmap=_cmap, vmin=-d_abs, vmax=d_abs,
                                   interpolation="nearest", aspect="equal")

                ax.axis("off")
                if row_idx == 0:
                    ax.set_title(title_txt, fontsize=8.5,
                                 fontweight="bold", pad=4)
                else:
                    ax.set_title(title_txt, fontsize=8, pad=4)
                panel_label(ax, LETTERS[letter_idx]); letter_idx += 1
                colorbar_h(fig, ax, im,
                           UNITS[var] if row_idx < 2 else f"Δ {UNITS[var]}")

        for row, rtitle in enumerate(row_titles):
            axes[row, 0].text(
                -0.08, 0.5, rtitle,
                transform=axes[row, 0].transAxes,
                fontsize=8, va="center", ha="right",
                rotation=90, style="italic",
            )

        fig.suptitle(
            f"Estilo Figura 6 — {var}: predicción GPR vs. relleno de gaps "
            f"en etapas fenológicas clave (Ciclo B, 2024)",
            fontsize=10, fontweight="bold", y=1.01,
        )
        out = os.path.join(OUTPUT_DIR, f"Fig6style_GPR_GF_{var}_CicloB.png")
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close()
        print(f"  ✅ {out}")


# ════════════════════════════════════════════════════════════
# FIGURA C: Serie temporal + barras fenológicas
# Etiquetas de etapa en eje superior PARA CICLO A Y CICLO B
# ════════════════════════════════════════════════════════════
def fig_timeseries_with_phenology():
    pheno_bg = [
        (0,   7,   0,   7,  "VE",     "#f7fcf0"),
        (7,  20,   7,  20,  "V1–V4",  "#e5f5e0"),
        (20, 35,  20,  35,  "V4–V8",  "#c7e9c0"),
        (35, 50,  35,  50,  "V8–V12", "#a1d99b"),
        (50, 65,  50,  65,  "V12–VT", "#74c476"),
        (65, 70,  65,  70,  "VT",     "#fdae6b"),
        (70, 75,  70,  75,  "R1",     "#fd8d3c"),
        (75, 85,  75,  85,  "R2",     "#f16913"),
        (85, 95,  85,  95,  "R3",     "#d94801"),
        (95,105,  95, 105,  "R4",     "#a63603"),
        (105,115,105, 115,  "R5",     "#7f2704"),
        (115,120,115, 120,  "R6",     "#4d1c00"),
    ]

    # DOY acumulados desde siembra para cada etiqueta del eje superior
    stage_offsets = [0, 7, 20, 35, 50, 65, 70, 75, 85, 95, 105, 115, 120]
    stage_lbls    = ["VE","V1","V4","V8","V12","VT","R1","R2","R3","R4","R5","R6",""]

    doy_a_start = 106   # 15 abr
    doy_b_start = 233   # 20 ago

    fig, axes = plt.subplots(
        len(VARIABLES), 1,
        figsize=(13, 4.0 * len(VARIABLES)),
        sharex=True,
        gridspec_kw={"hspace": 0.65}    # más espacio para dos ejes superiores
    )

    for ax_idx, (ax, var) in enumerate(zip(axes, VARIABLES)):
        folder_gf  = os.path.join(BASE_DIR, f"GF_{var}_2024")
        folder_gpr = os.path.join(BASE_DIR, f"GPR_{var}_2024")
        if not os.path.isdir(folder_gf): continue

        dates_gf,  arr_gf,  _ = load_collection(folder_gf)
        dates_gpr, arr_gpr, _ = load_collection(folder_gpr)

        med_gf  = [np.nanmean(np.where(a > 0, a, np.nan)) for a in arr_gf]
        p25_gf  = [np.nanpercentile(np.where(a>0,a,np.nan), 25)
                   for a in arr_gf]
        p75_gf  = [np.nanpercentile(np.where(a>0,a,np.nan), 75)
                   for a in arr_gf]
        med_gpr = [np.nanmean(np.where(a > 0, a, np.nan))
                   if np.any(a > 0) else np.nan for a in arr_gpr]

        doys_gf  = [d.timetuple().tm_yday for d in dates_gf]
        doys_gpr = [d.timetuple().tm_yday for d in dates_gpr]

        # Fondo fenológico
        for d0, d1, d0b, d1b, lbl, col in pheno_bg:
            ax.axvspan(doy_a_start + d0, doy_a_start + d1,
                       alpha=0.20, color=col, lw=0)
            ax.axvspan(doy_b_start + d0b, doy_b_start + d1b,
                       alpha=0.20, color=col, lw=0)

        ax.axvline(doy_a_start, color="#888", lw=0.8,
                   linestyle="--", alpha=0.6)
        ax.axvline(doy_b_start, color="#888", lw=0.8,
                   linestyle="--", alpha=0.6)

        # Curvas
        ax.fill_between(doys_gf, p25_gf, p75_gf,
                        alpha=0.18, color="#E07B3F")
        ax.plot(doys_gf, med_gf, color="#E07B3F", lw=1.8, zorder=3,
                label="GF relleno de gaps (mediana)")
        ax.scatter(doys_gpr, med_gpr, color="#2166AC", s=25,
                   edgecolors="white", lw=0.4, zorder=4,
                   label="Estimación GPR S2 (mediana)")

        ax.set_ylabel(f"{var}\n({UNITS[var]})", fontsize=9)
        ax.legend(loc="upper left", fontsize=7.5, framealpha=0.85)
        ax.grid(True, linestyle=":", lw=0.4, alpha=0.5)

        # ── Eje superior Ciclo A ──────────────────────────────
        ax2 = ax.twiny()
        ticks_a = [doy_a_start + d for d in stage_offsets]
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(ticks_a)
        ax2.set_xticklabels(stage_lbls, fontsize=6.5, color="#1a5276")
        ax2.tick_params(length=3, width=0.5, color="#1a5276",
                        labelcolor="#1a5276", top=True,
                        labeltop=True, direction="out")
        ax2.spines["top"].set_edgecolor("#1a5276")
        ax2.spines["top"].set_linewidth(0.6)

        # Etiqueta "Ciclo A" sobre el eje superior
        ax2.text(0.13, 1.045, "Ciclo A →",
                 transform=ax2.transAxes,
                 fontsize=7, color="#1a5276",
                 fontweight="bold", va="bottom", ha="left")

        # ── Eje superior Ciclo B ──────────────────────────────
        ax3 = ax.twiny()
        ticks_b = [doy_b_start + d for d in stage_offsets]
        ax3.set_xlim(ax.get_xlim())
        ax3.set_xticks(ticks_b)
        ax3.set_xticklabels(stage_lbls, fontsize=6.5, color="#7d6608")
        # Desplazar este eje un poco más arriba para que no se solape con ax2
        ax3.spines["top"].set_position(("axes", 1.10))
        ax3.tick_params(length=3, width=0.5, color="#7d6608",
                        labelcolor="#7d6608", top=True,
                        labeltop=True, direction="out")
        ax3.spines["top"].set_edgecolor("#7d6608")
        ax3.spines["top"].set_linewidth(0.6)

        # Etiqueta "Ciclo B" sobre el segundo eje superior
        ax3.text(0.60, 1.045, "Ciclo B →",
                 transform=ax3.transAxes,
                 fontsize=7, color="#7d6608",
                 fontweight="bold", va="bottom", ha="left")

        # Etiquetas de inicio de ciclo en el área del gráfico
        ymax_approx = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 1
        ax.text(doy_a_start + 2, ymax_approx * 0.96,
                "Ciclo A", fontsize=7, color="#1a5276",
                va="top", style="italic", fontweight="bold")
        ax.text(doy_b_start + 2, ymax_approx * 0.96,
                "Ciclo B", fontsize=7, color="#7d6608",
                va="top", style="italic", fontweight="bold")

    axes[-1].set_xlabel("Día del año (DOY) — 2024", fontsize=9)

    month_doys = {"Ene":1,"Feb":32,"Mar":60,"Abr":91,"May":121,
                  "Jun":152,"Jul":182,"Ago":213,"Sep":244,
                  "Oct":274,"Nov":305,"Dic":335}
    axes[-1].set_xticks(list(month_doys.values()))
    axes[-1].set_xticklabels(list(month_doys.keys()), fontsize=8)
    axes[-1].set_xlim(1, 366)

    fig.suptitle(
        "Perfiles temporales de LAI, FVC y laiCab — "
        "predicción GPR y relleno de gaps (2024)",
        fontsize=11, fontweight="bold", y=1.01
    )

    out = os.path.join(OUTPUT_DIR, "FigTS_timeseries_phenology.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  ✅ {out}")


# ════════════════════════════════════════════════════════════
# FIGURA D: LSP maps — panel 3×4
# POS y LOS: exterior blanco, texto leyenda negro
# ════════════════════════════════════════════════════════════
def fig_lsp_full_panel():
    LSP_BANDS  = ["sos", "eos", "pos", "los"]
    LSP_TITLES = {"sos": "SOS (DOY)", "eos": "EOS (DOY)",
                  "pos": "POS (DOY)", "los": "LOS (DOY)"}
    LSP_CMAPS  = {"sos": "YlOrBr", "eos": "YlOrBr",
                  "pos": "RdYlGn",  "los": "Blues"}
    CLIP_BANDS = {"pos", "los"}

    fig, axes = plt.subplots(
        len(VARIABLES), 4,
        figsize=(14, 3.8 * len(VARIABLES)),
        gridspec_kw={"hspace": 0.40, "wspace": 0.12}
    )
    # Fondo blanco para toda la figura (celdas vacías y exteriores al pol.)
    fig.patch.set_facecolor("white")

    letter_idx = 0

    for row, var in enumerate(VARIABLES):
        lsp_f = os.path.join(BASE_DIR, f"LSP_{var}_2024.tif")
        if not os.path.exists(lsp_f):
            for col in range(4): axes[row, col].set_visible(False)
            continue

        arrays      = {}
        transforms_ = {}
        crss_       = {}
        with rasterio.open(lsp_f) as src:
            nd      = src.nodata
            src_crs = src.crs
            src_tf  = src.transform
            for bi in range(src.count):
                bname = LSP_BANDS[bi] if bi < 4 else f"b{bi}"
                desc  = (src.descriptions[bi] or "").lower()
                key   = next((b for b in LSP_BANDS if b in desc), bname)
                arr   = clean(src.read(bi+1).astype(np.float32), nd)
                arrays[key]      = arr
                transforms_[key] = src_tf
                crss_[key]       = src_crs

        ref_arr = list(arrays.values())[0]
        extent  = get_aoi(ref_arr, pad=5)

        for col, bname in enumerate(LSP_BANDS):
            ax = axes[row, col]
            ax.set_facecolor("white")      # fondo del axes = blanco

            if bname not in arrays:
                ax.set_visible(False)
                letter_idx += 1
                continue

            arr_raw = arrays[bname]

            # ── Recorte estricto al polígono para POS y LOS ──
            if bname in CLIP_BANDS:
                arr_masked = mask_array_to_finca(
                    arr_raw, transforms_[bname], crss_[bname])
            else:
                arr_masked = arr_raw

            arr_c = crop(arr_masked, extent)

            # Usar masked array → NaN se renderiza con set_bad("white")
            arr_ma = np.ma.masked_invalid(arr_c)

            valid = arr_ma.compressed()
            vmin  = float(np.percentile(valid, 2))  if valid.size > 0 else 0
            vmax  = float(np.percentile(valid, 98)) if valid.size > 0 else 1

            # Copia del cmap con NaN → blanco
            import copy
            _cmap = copy.copy(plt.cm.get_cmap(LSP_CMAPS[bname]))
            _cmap.set_bad(color="white")    # ← exterior blanco
            _cmap.set_under(color="white")

            im = ax.imshow(arr_ma, cmap=_cmap,
                           vmin=vmin, vmax=vmax,
                           interpolation="nearest", aspect="equal")

            ax.axis("off")
            panel_label(ax, LETTERS[letter_idx]); letter_idx += 1

            # Colorbar con texto negro explícito
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("bottom", size="6%", pad=0.10)
            cb  = fig.colorbar(im, cax=cax, orientation="horizontal")
            cb.set_label("DOY" if bname != "los" else "días",
                         fontsize=7, labelpad=1, color="black")
            cb.ax.tick_params(labelsize=7, labelcolor="black",
                              color="black")
            cb.outline.set_edgecolor("black")
            cb.locator = mticker.MaxNLocator(nbins=5)
            cb.update_ticks()

            if row == 0:
                ax.set_title(LSP_TITLES[bname],
                             fontsize=9, fontweight="bold",
                             pad=5, color="black")
            if col == 0:
                ax.text(-0.10, 0.5, var,
                        transform=ax.transAxes,
                        fontsize=10, fontweight="bold",
                        va="center", ha="right",
                        rotation=90, color="black")

    fig.suptitle(
        "Métricas de fenología de la superficie terrestre (LSP) — "
        "LAI, FVC, laiCab (2024)",
        fontsize=11, fontweight="bold", color="black"
    )
    out = os.path.join(OUTPUT_DIR, "FigLSP_full_panel.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  ✅ {out}")


# ─── EJECUCIÓN PRINCIPAL ─────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  Figuras artículo — estilo IJRS/MDPI")
    print("=" * 60)

    print("\n[1] Mapa de ubicación (3 paneles)...")
    fig_location_map()

    print("\n[2] Comparación GPR→GF estilo Figura 6 (Ciclo B)...")
    fig_gpr_gf_comparison()

    print("\n[3] Series temporales con fenología (Ciclo A + B)...")
    fig_timeseries_with_phenology()

    print("\n[4] Panel completo LSP (3×4)...")
    fig_lsp_full_panel()

    print(f"\n✅ Todas las figuras en: {OUTPUT_DIR}")
    