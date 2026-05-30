// Import the LSP functions
var phenoFuncs = require('users/msalinero85/GPRPhenologyDemos:PhenologyFunctions');
var addDoy = phenoFuncs.addDoy;
var get_Double_Logistic_Params = phenoFuncs.get_Double_Logistic_Params

// Set the index
var vegIndex = 'laiCab'; 

// Pre-generated gapfilled collection
var gapfilledCol = ee.ImageCollection("users/msalinero85/30TUM_2019_Gapfilled_" + vegIndex).map(addDoy);

// Obtain the LSP metrics
var LSPmetrics = get_Double_Logistic_Params(gapfilledCol,'gapfilled',0.3);

// Visualization
var lsp_palette = ['FFFFFF','8EF82D','78C135','3FBA37','DAEB0D','EBB60D','EB740D','F20C0C','E817DF','20A1E5','1835D7','49538C'];

var lspVis = {
  palette: lsp_palette,
  min : 0,
  max : 365
}

var lsp = 'eos'

Map.centerObject(gapfilledCol.first(),9.6);
Map.addLayer(LSPmetrics.select(lsp), lspVis, vegIndex+'_'+lsp);

//Export 
/*
Export.image.toAsset({
  image: LSPmetrics.select(['sos', 'eos', 'pos', 'los']),
  description:  'PHENO_30TUM_2019_' + vegIndex,
  assetId: 'PHENO_30TUM_2019_' + vegIndex,
  scale: 20,
  region: gapfilledCol.first().geometry(),
  maxPixels: 1e13
});
*/
