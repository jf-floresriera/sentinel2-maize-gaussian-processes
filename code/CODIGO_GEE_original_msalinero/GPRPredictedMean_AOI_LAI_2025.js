
// =============================================================================
// SECCIÓN 1 · ÁREA DE INTERÉS (AOI) — Finca La Esperanza sin Lote Casa
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

var AOI = Area_Total.difference({right: Lote_Casa, maxError: 1});
var geometry = AOI;
var parcel = AOI;

var area_total_ha = Area_Total.area(1).divide(10000);
var area_casa_ha  = Lote_Casa.area(1).divide(10000);
var area_neta_ha  = AOI.area(1).divide(10000);
print('Área Total finca (ha):', area_total_ha);
print('Área Lote_Casa (ha):', area_casa_ha);
print('Área finca sin casa (ha):', area_neta_ha);

Map.centerObject(AOI, 14);
Map.addLayer(Area_Total, {color: '98ff00'}, 'Área Total');
Map.addLayer(Lote_Casa, {color: 'ff0000'}, 'Lote Casa (excluido)');
Map.addLayer(AOI, {color: '0080ff', opacity: 0.4}, 'AOI - Finca sin Casa');


// ============================================================================
// GPRPredictedMean ajustado a AOI + LAI + 2025
// ============================================================================

// Veg Index
var currentVegIndex = 'LAI';

// BOA Models
var importedModels = require('users/msalinero85/GPRPhenologyDemos:S2BOAModels');
var currentModel = importedModels.models[currentVegIndex];

// Visualization
var vis = require('users/msalinero85/GPRPhenologyDemos:visualization');
var currentVis = vis.visparams[currentVegIndex];

// FUNCTIONS
//Masking function for water and cloud
function maskS2cloud_and_water(image){
  
  var not_saturated = image.select('SCL').neq(1);
  var not_darl_area = image.select('SCL').neq(2);
  var not_cloud_shadows = image.select('SCL').neq(3);
  var not_water = image.select('SCL').neq(6);
  var not_cloud_low = image.select('SCL').neq(7);
  var not_cloud_medium = image.select('SCL').neq(8);
  var not_cloud_high = image.select('SCL').neq(9);
  var not_cirrus = image.select('SCL').neq(10);
  var not_ice = image.select('SCL').neq(11);
  
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(qa.bitwiseAnd(cirrusBitMask).eq(0)).and(not_water).and(not_cloud_low)
             .and(not_cloud_medium).and(not_cloud_high).and(not_cirrus).and(not_cloud_shadows).and(not_ice);
  return image.updateMask(mask).divide(currentModel.scaleFactor).copyProperties(qa).set('system:time_start', qa.get('system:time_start'));
}

// Veg index GPR
var veg_index_GPR = function(image_orig){
  
  //Create List of Bands of Dimension n (Xtrain[n,n])
  var XTrain_dim = currentModel.X_train.length().get([0]);
  var band_sequence   = ee.List.sequence(1, XTrain_dim).map(function(element){ return ee.String('B').cat(ee.String(element)).replace('[.]+[0-9]*$','')});

  //Create a list of band names for flattening operation 
  var im_norm_ell2D_hypell = image_orig.subtract(ee.Image(currentModel.mx)).divide(ee.Image(currentModel.sx)).multiply(ee.Image(currentModel.hyp_ell)).toArray().toArray(1); 
  var im_norm_ell2D = image_orig.subtract(ee.Image(currentModel.mx)).divide(ee.Image(currentModel.sx)).toArray().toArray(1); 
  var PtTPt  = im_norm_ell2D_hypell.matrixTranspose().matrixMultiply(im_norm_ell2D).arrayProject([0]).multiply(-0.5); //OK
  
  var PtTDX  = ee.Image(currentModel.X_train).matrixMultiply(im_norm_ell2D_hypell).arrayProject([0]).arrayFlatten([band_sequence]);
  var arg1   = PtTPt.exp().multiply(currentModel.hyp_sig);
  var k_star = PtTDX.subtract(ee.Image(currentModel.XDX_pre_calc).multiply(0.5)).exp().toArray();
  var mean_pred = k_star.arrayDotProduct(ee.Image(currentModel.alpha_coefficients).toArray()).multiply(arg1);
  mean_pred = mean_pred.toArray(1).arrayProject([0]).arrayFlatten([[currentModel.veg_index]]);
  mean_pred = mean_pred.add(currentModel.mean_model);
  mean_pred = mean_pred.where(mean_pred.lt(0),ee.Image(0.00001))
  image_orig= image_orig.addBands(mean_pred);
  return image_orig.select(currentModel.veg_index);
};


// Image retrieval
var Date_Start_str = '2025-01-01';
var Date_End_str = '2025-12-01';


var veg_index_collection  = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                            .filterBounds(geometry)
                            .filterDate(Date_Start_str,Date_End_str)
                            .map(maskS2cloud_and_water)
                            .select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'])
                            .map(veg_index_GPR);


var pm_image = veg_index_collection.first();
Map.centerObject(geometry, 14);
Map.addLayer(pm_image,currentVis,currentVegIndex);


// Color bar 
function makeColorBarParams(palette) {
  return {
    bbox: [0, 0, 1, 0.1],
    dimensions: '300x20',
    format: 'png',
    min: 0,
    max: 1,
    palette: currentVis.palette,
  };
}

// Create the color bar for the legend.
var colorBar = ui.Thumbnail({
  image: ee.Image.pixelLonLat().select(0),
  params: makeColorBarParams(vis.palette),
  style: {stretch: 'horizontal', margin: '0px 8px', maxHeight: '24px'},
});

// Create a panel with three numbers for the legend.
var legendLabels = ui.Panel({
  widgets: [
    ui.Label(0, {margin: '4px 8px'}),
    ui.Label(
        (400),
        {margin: '4px 8px', textAlign: 'center', stretch: 'horizontal'}),
    ui.Label(
        (800),
        {margin: '4px 8px', textAlign: 'center', stretch: 'horizontal'}),
    ui.Label(
        (1200),
        {margin: '4px 8px', textAlign: 'center', stretch: 'horizontal'}),
    ui.Label(
        (1600),
        {margin: '4px 8px', textAlign: 'center', stretch: 'horizontal'}),
    ui.Label(2000, {margin: '4px 8px'})
  ],
  layout: ui.Panel.Layout.flow('horizontal')
});

var legendTitle = ui.Label({
  value: currentVegIndex + '(' +  currentModel.units + ')',
  style: {fontWeight: 'bold', textAlign: 'center', stretch: 'horizontal'}
});

// Add the legendPanel to the map.
var legendPanel = ui.Panel([legendTitle, colorBar, legendLabels]);
Map.add(legendPanel);


var RGB_image  = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                 .filterBounds(geometry)
                 .filterDate(Date_Start_str,Date_End_str)
                 .first()

var RGBvis = {
  min: 0,
  max: 3000,
  bands:['B4','B3','B2'],
};          

Map.addLayer(RGB_image,RGBvis,'RGB');


