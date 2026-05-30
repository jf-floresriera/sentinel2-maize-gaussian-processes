// ============================================================
// SCRIPT COMPLETO FINAL: GPR LAI / FVC / laiCab
// Finca La Esperanza — Año 2024
// SIN .getInfo() en loops → no congela el navegador
// ============================================================

var importedModels = require('users/msalinero85/GPRPhenologyDemos:S2BOAModels');
var vis            = require('users/msalinero85/GPRPhenologyDemos:visualization');

// =============================================================================
// SECCIÓN 1 · AOI
// =============================================================================
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

var AOI       = Area_Total.difference({ right: Lote_Casa, maxError: 1 });
var parcel    = AOI.bounds();
var maiz_mask = ee.Image(1).clip(AOI).rename('maiz');

print('Área Total finca (ha):',    Area_Total.area(1).divide(10000));
print('Área Lote_Casa (ha):',      Lote_Casa.area(1).divide(10000));
print('Área finca sin casa (ha):', AOI.area(1).divide(10000));

Map.centerObject(AOI, 14);
Map.addLayer(Area_Total, {color: '98ff00'},              'Área Total');
Map.addLayer(Lote_Casa,  {color: 'ff0000'},              'Lote Casa (excluido)');
Map.addLayer(AOI,        {color: '0080ff', opacity: 0.4},'AOI — Finca sin Casa');

// =============================================================================
// SECCIÓN 2 · CONFIGURACIÓN
// =============================================================================
var currentTile  = '18NZK';
var start_date   = '2024-01-01';
var end_date     = '2024-12-31';
var ASSET_FOLDER = 'projects/wide-origin-466923-d8/assets';
var vegIndices   = ['LAI', 'FVC', 'laiCab'];

// =============================================================================
// SECCIÓN 3 · MÁSCARA DE NUBES
// ⚠️ Solo .set() — NO copyProperties (devuelve Element y rompe la cadena)
// =============================================================================
function maskS2cloud_and_water(image) {
  var scl  = image.select('SCL');
  var qa   = image.select('QA60');
  var mask = scl.neq(1).and(scl.neq(2)).and(scl.neq(3))
               .and(scl.neq(6)).and(scl.neq(7)).and(scl.neq(8))
               .and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11))
               .and(qa.bitwiseAnd(1 << 10).eq(0))
               .and(qa.bitwiseAnd(1 << 11).eq(0));
  return image.updateMask(mask)
    .set('system:time_start', image.get('system:time_start'));
}

// =============================================================================
// SECCIÓN 4 · MODELO GPR
// =============================================================================
function makeGPR(model) {
  return function(image_orig) {
    var XTrain_dim = model.X_train.length().get([0]);
    var band_seq   = ee.List.sequence(1, XTrain_dim).map(function(e) {
      return ee.String('B').cat(ee.String(e)).replace('[.]+[0-9]*$', '');
    });
    var im_hyp  = image_orig
      .subtract(ee.Image(model.mx)).divide(ee.Image(model.sx))
      .multiply(ee.Image(model.hyp_ell)).toArray().toArray(1);
    var im_norm = image_orig
      .subtract(ee.Image(model.mx)).divide(ee.Image(model.sx))
      .toArray().toArray(1);
    var PtTPt  = im_hyp.matrixTranspose().matrixMultiply(im_norm)
                   .arrayProject([0]).multiply(-0.5);
    var PtTDX  = ee.Image(model.X_train)
                   .matrixMultiply(im_hyp)
                   .arrayProject([0]).arrayFlatten([band_seq]);
    var arg1   = PtTPt.exp().multiply(model.hyp_sig);
    var k_star = PtTDX.subtract(
                   ee.Image(model.XDX_pre_calc).multiply(0.5)).exp().toArray();
    var mean_pred = k_star
      .arrayDotProduct(ee.Image(model.alpha_coefficients).toArray())
      .multiply(arg1)
      .toArray(1).arrayProject([0])
      .arrayFlatten([[model.veg_index]])
      .add(model.mean_model);
    mean_pred = mean_pred.where(mean_pred.lt(0), ee.Image(0.00001));
    return image_orig.addBands(mean_pred)
      .select(model.veg_index).toFloat();
  };
}

// =============================================================================
// SECCIÓN 5 · addVariables
// =============================================================================
var addVariables = function(image) {
  var days = ee.Date(image.get('system:time_start'))
               .difference(ee.Date('1970-01-01'), 'days');
  return image.addBands(ee.Image(days).rename('t').float())
    .set('system:time_start', image.get('system:time_start'));
};

// =============================================================================
// SECCIÓN 6 · GAP-FILLING
// =============================================================================
function makeGapFill(model, ell2_ts, sigf_ts, sign_ts, veg_index_GPR_fn) {
  return function(image) {
    var imageDate = ee.Date(image.get('system:time_start'));
    var tnum      = imageDate.difference(ee.Date('1970-01-01'), 'days');
    var t0        = imageDate.advance(-30, 'day');
    var t1        = imageDate.advance( 30, 'day');

    var col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
      .filter(ee.Filter.eq('MGRS_TILE', currentTile))
      .filterBounds(parcel)
      .filterDate(t0, t1)
      .map(function(img) {
        var m = maskS2cloud_and_water(img);
        var s = m.divide(model.scaleFactor);
        var c = s.clip(parcel);
        return c.set('system:time_start', img.get('system:time_start'));
      })
      .select(['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'])
      .map(veg_index_GPR_fn)
      .map(function(img) {
        return img.updateMask(maiz_mask)
          .set('system:time_start', img.get('system:time_start'));
      });

    var mask_col = col.map(function(img) {
      return img.select(model.veg_index).gt(0).rename('msk')
        .set('system:time_start', img.get('system:time_start'));
    });
    var t_col   = col.map(addVariables).select('t');
    var t_star  = t_col.first().multiply(0).add(tnum).toArray().toArray(1);
    var N       = t_col.size();
    var t_vec   = mask_col.toBands().multiply(t_col.toBands())
                    .unmask().toArray().toArray(1);
    var v_vec   = col.toBands().unmask().toArray().toArray(1);
    var ones_v  = t_vec.multiply(0).add(1.0);
    var ones_s  = t_star.multiply(0).add(1.0);
    var II      = ee.Image(ee.Array.identity(N));
    var prod    = t_vec.matrixMultiply(ones_v.matrixTranspose());
    var K       = prod.subtract(prod.matrixTranspose())
                    .pow(2).multiply(ell2_ts).multiply(-0.5).exp().multiply(sigf_ts);
    var L       = II.multiply(sign_ts).add(K).matrixCholeskyDecomposition();
    var alpha   = L.matrixTranspose().matrixInverse()
                    .matrixMultiply(L.matrixInverse().matrixMultiply(v_vec));
    var t_s_mat = t_star.matrixMultiply(ones_v.matrixTranspose());
    var t_t_mat = t_vec.matrixMultiply(ones_s.matrixTranspose()).matrixTranspose();
    var k_s     = t_s_mat.subtract(t_t_mat)
                    .pow(2).multiply(ell2_ts).multiply(-0.5).exp().multiply(sigf_ts);
    var pred    = k_s.matrixMultiply(alpha).arrayProject([0])
                    .arrayFlatten([['gapfilled']]).toFloat();
    pred = pred.where(pred.lt(0), ee.Image(0.00001)).toFloat();
    return image.addBands(pred).select('gapfilled')
      .set('system:time_start', image.get('system:time_start'));
  };
}

// =============================================================================
// SECCIÓN 7 · COLECCIÓN BASE
// =============================================================================
var S2base = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filter(ee.Filter.eq('MGRS_TILE', currentTile))
  .filterBounds(parcel)
  .filterDate(start_date, end_date)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80));

print('Total imágenes S2 base 2024:', S2base.size());

Map.addLayer(
  S2base.filterDate('2024-03-01','2024-05-01')
    .map(maskS2cloud_and_water).median(),
  {bands: ['B4','B3','B2'], min: 0, max: 3000},
  'RGB abr-may 2024'
);

// =============================================================================
// SECCIÓN 8 · LOOP PRINCIPAL — genera exports con .evaluate() (NO .getInfo())
// ✅ evaluate() es ASÍNCRONO → no bloquea ni congela el navegador
// =============================================================================

// Función que genera los exports para UNA variable dado su índice numérico
function processVariable(vegIdx, model, visParams, ell2, sigf, sign) {

  var veg_index_GPR = makeGPR(model);

  var S2col = S2base
    .map(function(img) {
      var m = maskS2cloud_and_water(img);
      var s = m.divide(model.scaleFactor);
      var c = s.clip(parcel);
      return c.set('system:time_start', img.get('system:time_start'));
    })
    .select(['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'])
    .map(veg_index_GPR)
    .map(function(img) {
      return img.updateMask(maiz_mask)
        .set('system:time_start', img.get('system:time_start'));
    });

  print('── ' + vegIdx + ' · Imágenes 2024:', S2col.size());
  Map.addLayer(S2col.median(), visParams, vegIdx + ' — mediana 2024');

  // ── PASO 1: Export sin gap-filling ──────────────────────────────────
  // .evaluate() llama al servidor UNA SOLA VEZ y recibe todas las fechas
  // en el callback — sin bloquear el navegador
  var folder1 = ASSET_FOLDER + '/GPR_' + vegIdx + '_2024';
  var dates1  = S2col.aggregate_array('system:time_start')
                  .map(function(t){ return ee.Date(t).format('YYYY-MM-dd'); });
  var list1   = S2col.toList(S2col.size());

  dates1.evaluate(function(dateArr) {
    for (var i = 0; i < dateArr.length; i++) {
      Export.image.toAsset({
        image:       ee.Image(list1.get(i)),
        description: vegIdx + '_' + dateArr[i],
        assetId:     folder1 + '/' + vegIdx + '_' + dateArr[i],
        region:      parcel,
        scale:       20,
        crs:         'EPSG:4326',
        maxPixels:   1e9
      });
    }
    print('✅ ' + vegIdx + ' — ' + dateArr.length + ' tasks sin GF en panel Tasks');
  });

  // ── PASO 2: Gap-filling + Export ─────────────────────────────────────
  var folder2      = ASSET_FOLDER + '/GF_' + vegIdx + '_2024';
  var gapFillFn    = makeGapFill(model, ell2, sigf, sign, veg_index_GPR);
  var gapfilledCol = S2col.map(gapFillFn);
  var list2        = gapfilledCol.toList(gapfilledCol.size());
  var dates2       = gapfilledCol.aggregate_array('system:time_start')
                      .map(function(t){ return ee.Date(t).format('YYYYMMdd'); });

  dates2.evaluate(function(dateArr) {
    for (var j = 0; j < dateArr.length; j++) {
      Export.image.toAsset({
        image:       ee.Image(list2.get(j)),
        description: 'GF_' + vegIdx + '_Meta_' + dateArr[j],
        assetId:     folder2 + '/' + dateArr[j],
        scale:       20,
        region:      parcel,
        maxPixels:   1e13
      });
    }
    print('✅ ' + vegIdx + ' — ' + dateArr.length + ' tasks con GF en panel Tasks');
  });
}

// Llamar la función para cada variable con IIFE para cerrar correctamente
for (var v = 0; v < vegIndices.length; v++) {
  (function(idx) {
    var vegIdx  = vegIndices[idx];
    var model   = importedModels.models[vegIdx];
    var visP    = vis.visparams[vegIdx];
    var ell2    = model.gf_hyperparams['media'].ell2_ts;
    var sigf    = model.gf_hyperparams['media'].sigf_ts;
    var sign    = model.gf_hyperparams['media'].sign_ts;
    processVariable(vegIdx, model, visP, ell2, sigf, sign);
  })(v);
}

print('⏳ Generando tasks de forma asíncrona...');
print('   Espera unos segundos y aparecerán en el panel Tasks →');
print('   Luego haz clic en RUN ALL');
