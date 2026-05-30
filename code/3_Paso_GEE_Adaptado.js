// ============================================================
// PASO 3 CORREGIDO: lee assets desde Folder (IndexedFolder)
// usando ee.data.listAssets() en lugar de ee.ImageCollection()
// ============================================================

var phenoFuncs             = require('users/msalinero85/GPRPhenologyDemos:PhenologyFunctions');
var addDoy                 = phenoFuncs.addDoy;
var get_Double_Logistic_Params = phenoFuncs.get_Double_Logistic_Params;

var ASSET_FOLDER = 'projects/wide-origin-466923-d8/assets';
var vegIndices   = ['LAI', 'FVC', 'laiCab'];

var lsp_palette = [
  'FFFFFF','8EF82D','78C135','3FBA37',
  'DAEB0D','EBB60D','EB740D','F20C0C',
  'E817DF','20A1E5','1835D7','49538C'
];
var lspVis = { palette: lsp_palette, min: 0, max: 365 };

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

Map.centerObject(AOI, 14);
Map.addLayer(AOI, {color: '0080ff', opacity: 0.4}, 'AOI — Finca sin Casa');

// ── Función clave: construye ImageCollection desde un Folder ──────────────
// Lee los assets del folder con listAssets(), los convierte a ee.Image
// y los empaqueta en una ImageCollection usando ee.ImageCollection(list)
function folderToImageCollection(folderPath) {
  // listAssets devuelve un objeto JS con campo 'assets' (array de {id, type})
  var assets   = ee.data.listAssets(folderPath);
  var assetIds = assets['assets'].map(function(a) { return a['id']; });
  var imgList  = assetIds.map(function(id) { return ee.Image(id); });
  return ee.ImageCollection(imgList);
}

// ── Loop por variable ─────────────────────────────────────────────────────
for (var v = 0; v < vegIndices.length; v++) {
  (function(vegIdx) {

    var folderPath = ASSET_FOLDER + '/GF_' + vegIdx + '_2024';

    // Construye la colección desde el Folder (no desde ImageCollection)
    var gapfilledCol = folderToImageCollection(folderPath).map(addDoy);

    print('── ' + vegIdx + ' · Imágenes gap-filled cargadas:', gapfilledCol.size());

    // Calcula métricas LSP
    var LSPmetrics = get_Double_Logistic_Params(gapfilledCol, 'gapfilled', 0.3);

    // Visualiza
    ['sos', 'eos', 'pos', 'los'].forEach(function(lsp) {
      Map.addLayer(
        LSPmetrics.select(lsp).clip(AOI),
        lspVis,
        vegIdx + ' — ' + lsp.toUpperCase() + ' 2024'
      );
    });

    print('Métricas LSP para ' + vegIdx + ':', LSPmetrics.bandNames());

    // Export
    Export.image.toAsset({
      image:       LSPmetrics.select(['sos', 'eos', 'pos', 'los']).clip(AOI),
      description: 'LSP_' + vegIdx + '_FincaLaEsperanza_2024',
      assetId:     ASSET_FOLDER + '/LSP_' + vegIdx + '_2024',
      scale:       20,
      region:      parcel,
      crs:         'EPSG:4326',
      maxPixels:   1e13
    });

    print('✅ ' + vegIdx + ' — Task LSP lista en panel Tasks');

  })(vegIndices[v]);
}

print('👉 Ve a Tasks → RUN ALL');
