var importedModels = require('users/msalinero85/GPRPhenologyDemos:S2BOAModels');

var vis = require('users/msalinero85/GPRPhenologyDemos:visualization');

var currentVegIndex = 'LAI';

var currentModel = importedModels.models[currentVegIndex];
var currentVis = vis.visparams[currentVegIndex];


// GPR hyperparameters for time series gapfilling
var ell2_ts  =  currentModel.gf_hyperparams['media'].ell2_ts;
var sigf_ts  =  currentModel.gf_hyperparams['media'].sigf_ts;
var sign_ts  =  currentModel.gf_hyperparams['media'].sign_ts;

/// FUNCTIONS

//Add a band with the days from epoch
var addVariables = function(image){
  var date = ee.Date(image.get("system:time_start"));
  var years = date.difference(ee.Date('1970-01-01'),'days');
  return image.addBands(ee.Image(years).rename('t').float());
};

//Masking function for water and cloud
function maskS2cloud_and_water(image){
  var not_water = image.select('SCL').neq(6);
  var not_cloud_low = image.select('SCL').neq(7);
  var not_cloud_medium = image.select('SCL').neq(8);
  var not_cloud_high = image.select('SCL').neq(9);
  var not_cirrus = image.select('SCL').neq(10);
  
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(qa.bitwiseAnd(cirrusBitMask).eq(0)).and(not_water).and(not_cloud_low)
             .and(not_cloud_medium).and(not_cloud_high).and(not_cirrus);
  return image.updateMask(mask).divide(currentModel.scaleFactor).copyProperties(qa).set('system:time_start', qa.get('system:time_start'));
}

var veg_index_GPR = function(image_orig){
  
  //Create List of Bands of Dimension n (Xtrain[n,n])
  var XTrain_dim = currentModel.X_train.length().get([0]);
  var band_sequence   = ee.List.sequence(1, XTrain_dim).map(function(element){ return ee.String('B').cat(ee.String(element)).replace('[.]+[0-9]*$','')});

  //Create a list of band names for flattening operation 
  var im_norm_ell2D_hypell = image_orig.subtract(ee.Image(currentModel.mx)).divide(ee.Image(currentModel.sx)).multiply(ee.Image(currentModel.hyp_ell)).toArray().toArray(1); 
  var im_norm_ell2D = image_orig.subtract(ee.Image(currentModel.mx)).divide(ee.Image(currentModel.sx)).toArray().toArray(1); 
  var PtTPt  = im_norm_ell2D_hypell.matrixTranspose().matrixMultiply(im_norm_ell2D).arrayProject([0]).multiply(-0.5); //OK
  
  //var PtTDX  = ee.Image(ee.Array(model.X_train)).matrixMultiply(im_norm_ell2D_hypell).arrayProject([0]).arrayFlatten([band_sequence]);
  var PtTDX  = ee.Image(currentModel.X_train).matrixMultiply(im_norm_ell2D_hypell).arrayProject([0]).arrayFlatten([band_sequence]);
  var arg1   = PtTPt.exp().multiply(currentModel.hyp_sig);
  var k_star = PtTDX.subtract(ee.Image(currentModel.XDX_pre_calc).multiply(0.5)).exp().toArray();
  var mean_pred = k_star.arrayDotProduct(ee.Image(currentModel.alpha_coefficients).toArray()).multiply(arg1);
  mean_pred = mean_pred.toArray(1).arrayProject([0]).arrayFlatten([[currentModel.veg_index]]);
  mean_pred = mean_pred.add(currentModel.mean_model);
  mean_pred = mean_pred.where(mean_pred.lt(0),ee.Image(0.00001))
  image_orig= image_orig.addBands(mean_pred);
  return image_orig.select(currentModel.veg_index).toFloat();
};


//
var calculate_Tile_gapfilling = function(image){
   
  var imageDate = ee.Date(image.get('system:time_start'));
  
  var time_span_days      = 30;
  
  var dateOfInterest      = imageDate.format('YYYY-MM-dd');
  var Date_Start_str      = ee.Date(dateOfInterest).advance(-time_span_days,'day');
  var Date_End_str        = ee.Date(dateOfInterest).advance(time_span_days,'day');
  var dateOfInterest_num  = ee.Date(dateOfInterest).difference(ee.Date('1970-01-01'),'days');
  
  var veg_index_collection  = ee.ImageCollection('COPERNICUS/S2_SR')
                              .select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12','QA60','SCL'])
                              .filterBounds(parcel)
                              .filterDate(Date_Start_str,Date_End_str)
                              .map(maskS2cloud_and_water)
                              .select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'])
                              .map(veg_index_GPR);
                                             
  var veg_index_collection_ini_msk =  veg_index_collection.map(function(image){ var veg_index_mask = image.select(currentModel.veg_index).gt(0).rename('veg_index_mask');  return image.addBands(veg_index_mask).select('veg_index_mask')});
  
  var time_im_Im_ini          =  veg_index_collection.map(addVariables).select('t');

  var t_star_vec              =  time_im_Im_ini.first().multiply(0).add(dateOfInterest_num).toArray().toArray(1);
  
  var Nsize       = time_im_Im_ini.size();
  
  var t_vec_sel  = veg_index_collection_ini_msk.toBands().multiply(time_im_Im_ini.toBands()).unmask().toArray().toArray(1);

  // LAI time series converted to vector+
  var lai_vec_sel     = veg_index_collection.toBands().unmask().toArray().toArray(1);
  
  var ones_vec_sel  = t_vec_sel.multiply(0).add(1.0);
  var ones_vec_star = t_star_vec.multiply(0).add(1.0);
  
  var II_mat        = ee.Image(ee.Array.identity(Nsize)); 
  var prod          = t_vec_sel.matrixMultiply(ones_vec_sel.matrixTranspose());
  var K_mat         = prod.subtract(prod.matrixTranspose()).pow(2).multiply(ell2_ts).multiply(-0.5).exp().multiply(sigf_ts);

  var L_mat        = II_mat.multiply(sign_ts).add(K_mat).matrixCholeskyDecomposition();
  var alpha_vec_1  = L_mat.matrixInverse().matrixMultiply(lai_vec_sel);
  var alpha_vec    = L_mat.matrixTranspose().matrixInverse().matrixMultiply(alpha_vec_1);
 
  var t_star_mat    = t_star_vec.matrixMultiply(ones_vec_sel.matrixTranspose());
  var t_train_mat   = t_vec_sel.matrixMultiply(ones_vec_star.matrixTranspose()).matrixTranspose();
  var k_star_mat    = t_star_mat.subtract(t_train_mat).pow(2).multiply(ell2_ts).multiply(-0.5).exp().multiply(sigf_ts);
  var pred_vec_aux  = k_star_mat.matrixMultiply(alpha_vec).arrayProject([0]);
  var pred_vec      = pred_vec_aux.arrayFlatten([['gapfilled']]).toFloat();
  pred_vec = pred_vec.where(pred_vec.lt(0),ee.Image(0.00001)).toFloat();
  image = image.addBands(pred_vec);
  return image.select('gapfilled');
  
};


// Prepare the datasets
var currentTile = '30TUM'
var start_date = '2019-04-11'
var end_date = '2019-04-15'

// Raw masked dataset
var S2collection_ini = ee.ImageCollection('COPERNICUS/S2_SR')
                       .filterBounds(parcel)
                       .filterDate(start_date,end_date)
                       .map(maskS2cloud_and_water)
                       .select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']);
                         
  
var img_gf = ee.Image(calculate_Tile_gapfilling(S2collection_ini.first()))
                       
// To avoid memory issues, it is highly recommended to export the image(as an asset or .tif) 
// and then visualize it. 
Map.centerObject(img_gf, 7.2);
Map.addLayer(img_gf,currentVis,currentVegIndex);


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


