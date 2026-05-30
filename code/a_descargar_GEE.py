# ============================================================
# Script a: Descarga de datos GEE imágenes descargadas
# ============================================================

import ee
import geemap
import os

# ─── AUTENTICACIÓN ───────────────────────────────────────────────────────────
ee.Initialize(project='wide-origin-466923-d8')

BASE_DIR = os.path.expanduser("~/GEE_Downloads")
os.makedirs(BASE_DIR, exist_ok=True)

PROJECT = "projects/wide-origin-466923-d8/assets"

IMAGE_COLLECTIONS = [
    (f"{PROJECT}/GPR_LAI_2024",    "GPR_LAI_2024"),
    (f"{PROJECT}/GF_LAI_2024",     "GF_LAI_2024"),
    (f"{PROJECT}/GPR_FVC_2024",    "GPR_FVC_2024"),
    (f"{PROJECT}/GF_FVC_2024",     "GF_FVC_2024"),
    (f"{PROJECT}/GPR_laiCab_2024", "GPR_laiCab_2024"),
    (f"{PROJECT}/GF_laiCab_2024",  "GF_laiCab_2024"),
]

SINGLE_IMAGES = [
    (f"{PROJECT}/LSP_LAI_2024",    "LSP_LAI_2024.tif"),
    (f"{PROJECT}/LSP_FVC_2024",    "LSP_FVC_2024.tif"),
    (f"{PROJECT}/LSP_laiCab_2024", "LSP_laiCab_2024.tif"),
]

# ─── FUNCIÓN: Escala y CRS nativos ───────────────────────────────────────────
def get_native_projection(image: ee.Image):
    proj  = image.select(0).projection()
    scale = proj.nominalScale().getInfo()
    crs   = proj.crs().getInfo()
    return scale, crs

# ─── FUNCIÓN CLAVE: Listar assets dentro de una carpeta/IndexedFolder ────────
def list_images_in_folder(folder_asset_path: str) -> list[dict]:
    """
    Usa ee.data.listAssets() para obtener todos los assets (imágenes)
    dentro de una carpeta, independientemente de si es ImageCollection
    o IndexedFolder.
    Maneja paginación automática con nextPageToken.
    """
    assets = []
    page_token = None

    while True:
        params = {"parent": folder_asset_path}
        if page_token:
            params["pageToken"] = page_token

        response = ee.data.listAssets(params)
        assets.extend(response.get("assets", []))

        page_token = response.get("nextPageToken")
        if not page_token:
            break

    return assets

# ─── DESCARGA DE CARPETAS (IndexedFolder o ImageCollection) ──────────────────
def download_collection(asset_path: str, folder_name: str):
    out_dir = os.path.join(BASE_DIR, folder_name)
    os.makedirs(out_dir, exist_ok=True)

    # ← CAMBIO CLAVE: listar con listAssets en vez de ee.ImageCollection
    asset_list = list_images_in_folder(asset_path)
    n = len(asset_list)
    print(f"\n📦 {folder_name}  →  {n} imágenes encontradas")

    if n == 0:
        print(f"  ⚠️  Carpeta vacía o sin acceso: {asset_path}")
        return

    for i, asset in enumerate(asset_list):
        img_id   = asset["name"].split("/")[-1]   # último segmento del path
        filename = os.path.join(out_dir, f"{img_id}.tif")

        if os.path.exists(filename):
            print(f"  ⏭  Ya existe: {img_id}.tif — omitiendo")
            continue

        # Cargar imagen por su path completo
        img = ee.Image(asset["name"])

        scale, crs = get_native_projection(img)
        region     = img.geometry().bounds()

        print(f"  ⬇  [{i+1}/{n}] {img_id}  (scale={scale:.1f}m, CRS={crs})")
        try:
            geemap.ee_export_image(
                img,
                filename=filename,
                scale=scale,
                crs=crs,
                region=region,
                file_per_band=False,
            )
        except Exception as e:
            print(f"  ❌ Error en {img_id}: {e}")

# ─── DESCARGA DE IMÁGENES LSP (multi-banda) ───────────────────────────────────
def download_single_image(asset_path: str, filename: str):
    out_path = os.path.join(BASE_DIR, filename)

    if os.path.exists(out_path):
        print(f"  ⏭  Ya existe: {filename} — omitiendo")
        return

    img        = ee.Image(asset_path)
    scale, crs = get_native_projection(img)
    region     = img.geometry().bounds()
    band_names = img.bandNames().getInfo()

    print(f"\n🗺  {filename}  →  bandas: {band_names}  (scale={scale:.1f}m, CRS={crs})")
    try:
        geemap.ee_export_image(
            img,
            filename=out_path,
            scale=scale,
            crs=crs,
            region=region,
            file_per_band=False,
        )
        print(f"  ✅ Descargado: {filename}")
    except Exception as e:
        print(f"  ❌ Error en {filename}: {e}")

# ─── MAIN ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("  GEE → Disco Local  |  Proyecto: wide-origin-466923-d8")
    print("=" * 60)

    for asset_path, folder_name in IMAGE_COLLECTIONS:
        download_collection(asset_path, folder_name)

    print("\n─── Imágenes LSP (multi-banda) ───")
    for asset_path, filename in SINGLE_IMAGES:
        download_single_image(asset_path, filename)

    print("\n✅ ¡Descarga completa!")
    print(f"   Archivos en: {BASE_DIR}")
    