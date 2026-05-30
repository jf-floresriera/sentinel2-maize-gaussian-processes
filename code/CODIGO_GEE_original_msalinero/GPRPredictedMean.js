// Veg Index
var currentVegIndex = 'laiCab'

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
var Date_Start_str = '2020-04-11'
var Date_End_str = '2020-04-15'


var veg_index_collection  = ee.ImageCollection('COPERNICUS/S2_SR')
                            .filterBounds(geometry)
                            .filterDate(Date_Start_str,Date_End_str)
                            .map(maskS2cloud_and_water)
                            .select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'])
                            .map(veg_index_GPR);


var pm_image = veg_index_collection.first();
Map.centerObject(pm_image,9.6);
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


var RGB_image  = ee.ImageCollection('COPERNICUS/S2_SR')
                 .filterBounds(geometry)
                 .filterDate(Date_Start_str,Date_End_str)
                 .first()

var RGBvis = {
  min: 0,
  max: 3000,
  bands:['B4','B3','B2'],
};          

Map.addLayer(RGB_image,RGBvis,'RGB');


