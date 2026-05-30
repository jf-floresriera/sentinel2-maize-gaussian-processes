//Get DOY of an image
exports.addDoy = function(image){
  var doy = image.date().getRelative('day', 'year');
  var doyBand = ee.Image.constant(doy).uint16().rename('doy');
  return image.addBands(doyBand).set('doy', doy);
};

//Generate LSP metrics
exports.get_Double_Logistic_Params = function(imageCol,targetBand,customGap){
  /*
  Function generating the necessary double logistic parameters.
  imageCol: imageCollection with the desired targetBand.
  targetBand: band from which extract the LPS metrics.
  customGap: custom threshold. 
  */
  
  //  0: get the vmin, vmax, and vamp
  var vmax = imageCol.filter(ee.Filter.and(ee.Filter.gte('doy', 60), ee.Filter.lte('doy', 304))).select(targetBand).reduce(ee.Reducer.intervalMean(95,100)).rename('vmax');
  var vmin = imageCol.filter(ee.Filter.and(ee.Filter.gte('doy', 60), ee.Filter.lte('doy', 304))).select(targetBand).reduce(ee.Reducer.intervalMean(0,5)).rename('vmin');  
  var vamp = vmax.subtract(vmin).rename('vamp');
  var vmid = vmax.add(vmin).multiply(0.5);
  
  /*  1: get gap2max and gap2mid (distance from the pixel to max and mid) and 
         obtain doymax (pixel with minimun distance to gaptomax)
  */
  imageCol = imageCol.map(function(image){
    var gap2max = image.select(targetBand).subtract(vmax).abs()
                  .rename('gap2max');
    var gap2mid = image.select(targetBand).subtract(vmid).abs()
                  .rename('gap2mid');
    return image.addBands([gap2max,gap2mid]);
  });
  var doymax = imageCol.select(['gap2max','doy']).reduce(ee.Reducer.min(2)).rename(['gap2max','doy']).select('doy');

  //  2: apply the doymax to findout two doymask for doyn1 and doyn2
  //     Limit doyn1 and doyn2 to be near doymax to avoid other seasons.
  var doyn_limit = 60
  imageCol = imageCol.map(function(image){
    var doy = image.select('doy');
  
    var valid11 = doy.subtract(doymax).lt(0);
    var valid12 = doy.subtract(doymax.subtract(doyn_limit)).gte(0);
    var valid1 = valid11.eq(1).add(valid12.eq(1)).eq(2);
    var gap2mid1 = image.select('gap2mid').updateMask(valid1).rename('gap2mid1');
    
    var valid21 = doy.subtract(doymax).gt(0);
    var valid22 = doy.subtract(doymax.add(doyn_limit)).lte(0);
    var valid2 = valid21.eq(1).add(valid22.eq(1)).eq(2);
    var gap2mid2 = image.select('gap2mid').updateMask(valid2).rename('gap2mid2');
    return image.addBands([gap2mid1,gap2mid2]);
  });
  
  var doyn1 = imageCol.select(['gap2mid1','doy']).reduce(ee.Reducer.min(2)).rename(['gap2mid','doy']).select('doy'); 
  var doyn2 = imageCol.select(['gap2mid2','doy']).reduce(ee.Reducer.min(2)).rename(['gap2mid','doy']).select('doy'); 
    
  //  3: using doyn1, doyn2 and doymax to update the imageCollection 
  var douLogFitCol1 = imageCol.map(function(image){
    var doymin = doyn1.multiply(ee.Image(2.0)).subtract(doymax);
    var valid1 = image.select('doy').subtract(doymin).gte(0); 
    var valid2 = image.select('doy').subtract(doymax).lte(0);
    var valid = valid1.eq(1).add(valid2.eq(1)).eq(2);
    return image.updateMask(valid);
  });
  
  var vmin1 = douLogFitCol1.select(targetBand).min().rename('vmin1');  
  var vamp1 = vmax.subtract(vmin1).rename('vamp1');
  
  var vsos = vamp1.multiply(customGap)
  vsos = vsos.add(vmin1);
  
  douLogFitCol1 = douLogFitCol1.map(function(image){
    var gap2sos = image.select(targetBand).subtract(vsos).abs()
                  .rename('gap2sos');
    return image.addBands(gap2sos);
  });
  
  var custom_sos = douLogFitCol1.select(['gap2sos','doy']).reduce(ee.Reducer.min(2)).rename(['gap2sos','doy']).select('doy').rename('custom_sos');

  var douLogFitCol2 = imageCol.map(function(image){
    
    var doymin = doyn2.multiply(ee.Image(2.0)).subtract(doymax);
    var valid1 = image.select('doy').subtract(doymin).lte(0); 
    var valid2 = image.select('doy').subtract(doymax).gte(0);
    var valid = valid1.eq(1).add(valid2.eq(1)).eq(2);
    return image.updateMask(valid);
    
  });
  
  
  var vmin2 = douLogFitCol2.select(targetBand).min().rename('vmin2');  
  var vamp2 = vmax.subtract(vmin2).rename('vamp2');
  
  var veos = vamp.multiply(customGap)
  veos = veos.add(vmin2);
  
  douLogFitCol2 = douLogFitCol2.map(function(image){
    var gap2eos = image.select(targetBand).subtract(veos).abs()
                  .rename('gap2eos');
    return image.addBands(gap2eos);
  });
  
  var custom_eos = douLogFitCol2.select(['gap2eos','doy']).reduce(ee.Reducer.min(2)).rename(['gap2eos','doy']).select('doy').rename('custom_eos');

  //  4: estimate the parameters 
  //  ****** (1) define function to get the n1, m1 (first sigmoid)
  douLogFitCol1 = douLogFitCol1.map(function(image){
    // convert the sigmoid function to linear form and get this parameter
    var sig_x = vmax.subtract(image.select(targetBand)).divide(image.select(targetBand).subtract(vmin)).log().rename('sig_x'); 
    var sig_y = image.select('doy').rename('sig_y');
    var constant = ee.Image.constant(1).float();
    
    return image.addBands([sig_x,sig_y,constant]);
  });
  
  //  apply linear regresssion for parameter estimation (n1, m1)
  var independents = ee.List(['constant', 'sig_x']);
  var dependent = ee.String('sig_y');
  var trend1 = ee.ImageCollection(douLogFitCol1).select(independents.add(dependent))
               .reduce(ee.Reducer.linearRegression(independents.length(), 1));
               
  var coefficients1 = trend1.select('coefficients').arrayProject([0]).arrayFlatten([independents]);
  var n1 = coefficients1.select(['constant']).float().rename('n1');
  var m1 = ee.Image(-1.0).divide(coefficients1.select(['sig_x'])).float().rename('m1');
    
  //  ****** (2) define function to get the n2, m2 (second sigmoid)
  douLogFitCol2 = douLogFitCol2.map(function(image){
    var tempX = image.select('doy');
    var tempSig1 = ee.Image(1.0).divide(ee.Image(1.0).add((m1.multiply(-1.0).multiply(tempX.subtract(n1))).exp()));
      
    // get the second sigmoid term
    var sig_x1 = vamp.multiply(ee.Image(1.0).subtract(tempSig1)).add(image.select(targetBand)).subtract(vmin);
    var sig_x2 = vamp.multiply(tempSig1).subtract(image.select(targetBand)).add(vmin);
    var sig_x = sig_x1.divide(sig_x2).log().rename('sig_x'); 
    var sig_y = image.select('doy').rename('sig_y');
    var constant = ee.Image.constant(1).float();
    return image.addBands([sig_x,sig_y,constant]);
  });
  var trend2 = ee.ImageCollection(douLogFitCol2).select(independents.add(dependent))
               .reduce(ee.Reducer.linearRegression(independents.length(), 1));
  var coefficients2 = trend2.select('coefficients').arrayProject([0]).arrayFlatten([independents]);
  var n2 = coefficients2.select(['constant']).float().rename('n2');
  var m2 = ee.Image(-1.0).divide(coefficients2.select(['sig_x'])).float().rename('m2');
  
  var pos = doymax.rename('pos');  
  var los = n2.subtract(n1).rename('los');
  var sos = n1.rename('sos');
  var eos = n2.rename('eos');
  
  var params = ee.Image().addBands([vmin, vmax, n1, m1, n2, m2, sos, eos, custom_sos, custom_eos, pos, los]).slice(1);//remove constant band
  return params;
};

exports.doubleLogistic_Fitting = function(DLParams, ic_doy){
  /*
  This function used the parameters in DLPara and the ImageCollection (ic_doy) 
  to get the prediction
  DLParams: n1, n2, m1, m2, vmin, vmax
  ic_doy: contain the DOY layer
  */
  
  //  get parameters
  var vmin = DLParams.select('vmin');
  var vmax = DLParams.select('vmax');
  var vamp = vmax.subtract(vmin); 
  var n1 = DLParams.select('n1');
  var n2 = DLParams.select('n2');
  var m1 = DLParams.select('m1');
  var m2 = DLParams.select('m2');
  
  return ee.ImageCollection(ic_doy).map(function(image){
    var tempX = image.select('doy');
    var term1 = tempX.subtract(n1).multiply(m1).multiply(ee.Image(-1.0)).exp().add(ee.Image(1)).pow(-1);
    var term2 = tempX.subtract(n2).multiply(m2).multiply(ee.Image(-1.0)).exp().add(ee.Image(1)).pow(-1);
    var tempY = term1.subtract(term2).multiply(vamp).add(vmin).rename('fitted').float(); 
    return image.addBands(tempY);
  })
}

