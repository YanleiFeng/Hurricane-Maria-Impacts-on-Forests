// Name: Yanlei Feng
// Last edit date: 01/12/2021
///////////////////////////////////////////////////////////////////
/////////////////////////***Functions***////////////////////////////
///////////////////////////////////////////////////////////////////
// Get BQA band bits
var getQABits = function(image, start, end, newName) {
    // Compute the bits we need to extract.
    var pattern = 0;
    for (var i = start; i <= end; i++) {
       pattern += Math.pow(2, i);
    }
    // Return a single band image of the extracted QA bits, giving the band
    // a new name.
    return image.select([0], [newName])
                  .bitwiseAnd(pattern)
                  .rightShift(start);
};

//make a buffer around clouds with a certain radius
var kernel = function(image){
  var kernelimage = image.focal_max({radius:1});
  return kernelimage;
} 

// A function to mask out cloudy pixels with kernel.
var maskClouds = function(image) {
  // Select the QA band.
  var QA = image.select('BQA');
  // Get the QA bits
  var internalCloud = getQABits(QA, 4, 4, 'L8PH_QA');
  // Find Cloudy pixels
  var cloud = internalCloud.eq(1)
  // Make a kernel around the pixels
  var cloudwkernel = kernel(cloud)
  var mask = cloudwkernel.unmask(0).clip(PRboundary).not();
  // Return an image masking out cloudy areas with kernels.
  var final = image.updateMask(mask);
  return final;
};

// A function to mask out cloud shadow pixels for post-hurricane images.
var maskCloudShadow = function(image) {
  // Select the QA band.
  var QA = image.select('BQA');
  // Get the QA bits
  var internalCloud = getQABits(QA, 7, 8, 'L8PH_CloudShadow');
  // Return an image masking out cloudy areas.
  return image.mask(internalCloud.eq(0).add(internalCloud.eq(1)));
};

//Mask all the other and leave only clear pixels in Landsat 8 Surface Reflectance
var clear = function(image){
  var QA = image.select('pixel_qa');
  var clearmask = (QA.eq(322).add(QA.eq(386)).add(QA.eq(834)).add(QA.eq(898)).add(QA.eq(1346))
                   .add(QA.eq(324)).add(QA.eq(388)).add(QA.eq(836)).add(QA.eq(900)).add(QA.eq(1348)))
  var clearimage = image.updateMask(clearmask)
  return clearimage;
}

//Mask Water for landsat 8
//Not in use
var maskWater = function(image) {
  var NDWI =  image.normalizedDifference(['B3','B5']).rename('NDWI')
  
  var watermask = NDWI.lte(0.5);
  var imagewithoutwater = image.mask(watermask);
  return imagewithoutwater.addBands(NDWI);
}

// Select Bands
function selectBands(image){
  return image.select(['B1','B2','B3','B4','B5','B6','B7']); // L8t1.
}

////////////////////////////////////////////////////
//Topographic Illumination by Jun Xiong, 06/20/2018
////////////////////////////////////////////////////
function illuminationCorr_namedBands (image){ 
var terrain = ee.call('Terrain', ee.Image('srtm90_v4'));
var solar_zenith = ee.Number((image.get("SOLAR_ZENITH_ANGLE")));
var solar_azimuth = ee.Number(image.get("SOLAR_AZIMUTH_ANGLE")); 
var degree2radian = 0.01745
var solar_zenith_radians = solar_zenith.multiply(degree2radian)
var slope_radians = terrain.select(['slope']).multiply(degree2radian)//expression("(b('slope')*0.01745");
var aspect = terrain.select(['aspect']);


//slope part of the illumination condition
var cosZ = solar_zenith_radians.cos();

var cosS = slope_radians.cos();
var slope_illumination = cosS.multiply(cosZ)//cosS.expression("b('slope')*(" + cosZ + ")").select(['slope'], ['b1']);

//aspect part of the illumination condition
var sinZ = solar_zenith_radians.sin();
var sinS = slope_radians.sin();
var azimuth_diff_radians = (aspect.subtract(solar_azimuth)).multiply(degree2radian)//.expression("((b('aspect')-" + solar_azimuth + ")*0.01745")

var cosPhi = azimuth_diff_radians.cos();
var aspect_illumination = cosPhi.multiply(sinS).multiply(sinZ)//.expression("b('aspect')*" + sinZ).select(['aspect'], ['b1']);

//illumination condition
var ic = slope_illumination.add(aspect_illumination)
//Map.addLayer(ic)
// Create the image to apply the linear regression.The first band
//is the ic and the second band is the response variable, the reflectance (the bands).
//L (y) = a + b*cosi(x); a = intercept, b = slope
// Dependent: Reflectance

var CoefficientLoop = function(band){
  //var y = image.select(band);
// Independent: (cosi)
  var x = ic.rename('ic');
// Intercept: a
  var a = ee.Image(1).rename('a');
// create an image collection with the three variables by concatenating them
  var total_img = image.addBands(x).addBands(a)
  var reg_img = total_img.select(['ic','a',band])//ee.Image.cat(x,a,y);
// specify the linear regression reducer
  var lr_reducer = ee.Reducer.linearRegression({
    numX: 2,
    numY: 1
});
// fit the model
  var fit = reg_img.reduceRegion({
   reducer: lr_reducer,
   geometry: TI,
   scale: 30,
   maxPixels: 1e10
});
//var fit = fit.combine({"coefficients": ee.Array([[1],[1]])}, false);

  var coef = ee.Array(fit.get("coefficients"))
  var slo = coef.get([0,0])
  var int = coef.get([1,0])
  var getC = int.divide(slo) //Calculate C parameter C= b/a, b is intercept and a is the slope
  var getA = slo
// Get the coefficients as a nested list, cast it to an array, and get
// just the selected column
//var slo = (ee.Array(fit.get('coefficients')).get([0,0]));
//var int = (ee.Array(fit.get('coefficients')).get([1,0]));


//var C = int.divide(slo);
  return (getC)
}

var bands = ['B1','B2','B3','B4','B5','B6','B7']
var Coefficients = bands.map(CoefficientLoop);
var A = bands.map(CoefficientLoop);



//apply the cosine correction
var cos_output = image.expression("((image*cosZ)/ic) + offsets", {
    'image': image.select(['B1','B2','B3','B4','B5','B6','B7']),
    'cosZ': cosZ,
    'ic': ic,
    'offsets': [0, 0, 0, 0, 0, 0, 0]
});

var c_output = image.expression("((image * (cosZ + coeff)) / (ic + coeff )) ", {
    'image': image.select(['B1','B2','B3','B4','B5','B6','B7']),
    'ic': ic,
    'cosZ': cosZ,
    'coeff': Coefficients,
    //'offsets': [0, 0, 0, 0, 0, 0, 0]
});

var c_output_scsc = image.expression(
    '((image * ((cosp*cosz) + C))/(ic + C))',
  {
      'image': image.select(['B1','B2','B3','B4','B5','B6','B7']),
      'cosp': slope_radians.cos(),
      'cosz': cosZ,
      'ic': ic,
      'C': Coefficients,
  });

var c_output_regress_tanEtAl = image.expression(
    'image - a * (ic - cosz)',
  {
      'image': image.select(['B1','B2','B3','B4','B5','B6','B7']),
      'cosz': cosZ,
      'ic': ic,
      'a': A,
  });

// print(c_output)
 return ee.Image( c_output_scsc
        .uint16()
        .addBands(image.select(['B10', 'B11','pixel_qa']))
        .copyProperties(image)
        .copyProperties(image, ['system:time_start']));
}



// Relative radiometric normalization, single band
// i_tbn: Image (single band) to be normalized
// i_b: Base image to use for normalization (single band)
// pir: pseudo invariant region (feature)
// returns: normalized image (single band)
var relNormalizeSB = function(i_tbn, i_b, pir) {

  var multiTemporal = i_tbn.addBands(ee.Image(1)).addBands(i_b)
  // Compute coefficients c1,c2 that minimizes error in c1*i_tbn + c2*1 = i_b
  var regression = multiTemporal.reduceRegion({
  reducer: ee.Reducer.robustLinearRegression(2, 1),
  geometry: pir, 
  scale: 30
  });
  // print(regression);
  // From dictionary, need to cast
  var coeffArray = ee.Array(regression.get('coefficients'))
  var c1 = coeffArray.get([0,0]);
  var c2 = coeffArray.get([1,0]);
  return (i_tbn.multiply(c1).add(c2));
};

// relative radiometric normalization, multiband
// i_tbn: Image (multiband) to be normalized.
// tbn: to be normalized
// i_b: Base image to use for normalization
// pir: pseudo invariant regions (feature)
// nb: number of bands to normalize. Eg, if nb=3,
// the first, second and third bands will be normalized
// and a three-band image will be returned.
// returns: normalized image
var relNormalize = function(i_tbn, i_b, pir, nb) {

  // the first band
  var rni = relNormalizeSB(i_tbn.select(0), i_b.select(0), pir);
  // Iterate over the rest of the bands
  for (var b = 1; b < nb; b++) {
    rni = rni.addBands(relNormalizeSB(i_tbn.select(b), i_b.select(b), pir));
  }
  return (rni).uint16()
}


var shade = function(image){
  var Endmember_shade = image.reduceRegions({
  collection: shade_water,
  reducer: ee.Reducer.mean(),
  scale: 30,
  });
  return Endmember_shade;
};


var gv = function(image){
  var Endmember_gv = image.reduceRegions({
  collection: shade_water,
  reducer: ee.Reducer.mean(),
  scale: 30,
  });
  return Endmember_shade;
};



// Spectral Mix Analysis with image-derived endmembers
var smaima =function(image){
var gv =   [285, 301, 771, 310, 7656, 2367, 853];   //GV: cecropia green @(-65.999476, 18.150572), Adam & Guillespie
var npv =  [375, 386, 500, 941, 2812, 4356, 2518];    //NPV
var shd =  [600, 763, 1595, 1877, 799, 394, 227]; //Shade: water Adam & Guillespie P273 @(-66.65543, 18.26408)
//var soil = [6119, 6828, 7884, 8567, 8317, 6273, 2570]//Soil (-66.55535, 18.02298). Add soil member result in very high DNPV
var sma = image.unmix([npv,gv,shd]);
return sma;
};

// Clip Function. clip the Puerto Rico boundary on images 
var clip = function(image) {
  return image.clip(PRboundary).clip(PRarea);
};


// // This function adds cloud quality bands to Landsat 8 images.
// var addQualityBands = function(image) {
//   return image
//     // Cloud Score
//     //.addBands((ee.Image(100).subtract(ee.Algorithms.Landsat.simpleCloudScore(image))), 'cloudscore')//simpleCloudScore only works for TOA
//     //NDVI
//     .addBands(image.normalizedDifference(['B5', 'B4']));
// };
// This function masks clouds and adds quality bands to Landsat 8 images.
var addQualityBands = function(image) {
  return image
    // NDVI
    .addBands(image.normalizedDifference(['B5', 'B4']))
    // time in days
    .addBands(image.metadata('system:time_start'));
};

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////***Main***/////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//Hurricane Maria made landfall in Puerto Rico on September 20. Category 4
var imacol=L8SR;
var filterarea = ee.Geometry.Rectangle(-67.3242, 17.9265, -65.517, 18.5421)
var PRcollection = imacol.filterBounds(filterarea)
                .filterDate("2017-10-01","2018-01-30")
                // .filterDate("2017-10-01","2017-11-01")
                // .filterDate("2018-06-01","2018-09-28")
               // .filterMetadata('CLOUD_COVER_LAND','less_than', 60).sort('CLOUD_COVER_LAND',false)
                .filterMetadata('CLOUD_COVER_LAND','less_than', 40).sort('CLOUD_COVER_LAND',false)
                
// var PRcollection = imacol.filterBounds(filterarea)
//                 .filterDate("2018-06-30","2018-09-30")
//                 //.filterDate("2017-10-01","2017-12-01")
//                 .filterMetadata('CLOUD_COVER_LAND','less_than', 50).sort('CLOUD_COVER_LAND',false)

var PRcollection2016 = imacol.filterBounds(filterarea)
                .filterDate("2016-06-01","2016-09-30") 
                .filterMetadata('CLOUD_COVER_LAND','less_than', 40).sort('CLOUD_COVER_LAND', false);

print("PRcollection2016", PRcollection2016);
//mask cloud and cloud shadow
var PostNoCloudnShadow = ee.ImageCollection(PRcollection.map(clear));
var PreNoCloudnShadow = ee.ImageCollection(PRcollection2016.map(clear));

Map.addLayer(PostNoCloudnShadow.map(clip), {bands: ['B6', 'B5', 'B4'], min:300, max: 5000},'Before TopoIllu: Post Hurricane')
// Map.addLayer(PreNoCloudnShadow.map(clip), {bands: ['B6', 'B5', 'B4'], min: 300, max: 5000},'Before TopoIllu: Pre Hurricane')

//Map Topographic Illumination on the ImageCollection
var PostNCSTI = PostNoCloudnShadow.map(illuminationCorr_namedBands);
var PreNCSTI = PreNoCloudnShadow.map(illuminationCorr_namedBands);
// Map.addLayer(PostNCSTI.map(clip), {bands: ['B6', 'B5', 'B4'], min:300, max: 5000},'After TopoIllu: Post Hurricane')
// Map.addLayer(PreNCSTI.map(clip), {bands: ['B6', 'B5', 'B4'], min: 300, max: 5000},'After TopoIllu: Pre Hurricane')

//Mosaic the least cloudy image on the top
var mosaic2017 = PostNCSTI.mosaic();
var mosaic2016 = PreNCSTI.mosaic();

// Display the composites.
Map.setCenter( -66.4838, 18.2033, 10);  
// Map.addLayer(mosaic2017, {bands: ['B6', 'B5', 'B4'], min: 300, max: 5000},
//     'mosaic2017');
// Map.addLayer(mosaic2016, {bands: ['B6', 'B5', 'B4'], min: 300, max: 5000},
//     'mosaic2016');

// image to be normalized (tbn: to be normalized)
var i_tbn = selectBands(mosaic2017);
// base image (to which the above image will be normalized to)
var i_base = selectBands(mosaic2016);
// normalize i_tbn image (first 7 bands only)
// Hand-selected 17 Invariant Targets, including 7 features in urban areas, 
//3 in forests, 4 in lakes, and 1 in farmland. 
//Invariant Targets were stored in the Geometry layer named 'IT', the targets are very small polygons
var i_norm = relNormalize(i_tbn, i_base, IT, 7);

// To compare images
Map.setCenter( -66.4838, 18.2033, 10);  
Map.addLayer(i_base, {bands: 'B6,B5,B4', 'gain':'0.06,0.04,0.06'}, 'Base: Pre-hurricane image');
// Map.addLayer(i_tbn, {bands: 'B6,B5,B4','gain':'0.06,0.04,0.06'}, 'Before RadioNorm: Post-hurricane image');
// Map.addLayer(i_norm, {bands: 'B6,B5,B4', 'gain':'0.06,0.04,0.06'}, 'After RadioNorm: Post-hurricane image.');


var shade = shade(i_base);
print('Shade endmembers', shade)
var studyarea = ee.Geometry.Rectangle(-67.35207, 17.87306, -65.54865, 18.56641)//lon: 303:345; lat: 308:325 in wind data
// To compare the mean reflectance value of bands in the invariant targets
// region. They should be almost same after normalization.
// print ('Base:');
// print(i_base.reduceRegion({
//   reducer: ee.Reducer.mean(),
//   geometry: IT,
//   scale:30
// }));
// Export.image.toDrive({
//   image: i_base,
//   description: 'PrehurricaneImage',
//   folder:'PRlndst8',
//   scale: 30,
//   region: studyarea
// });
// print ('Original image:')
// print(i_tbn.reduceRegion({
//   reducer: ee.Reducer.mean(),
//   geometry: IT,
//   scale:30
// }));
// print ('Normalized image:')
// print(i_norm.reduceRegion({
//   reducer: ee.Reducer.mean(),
//   geometry: IT,
//   scale:30
// }));

// Export.image.toDrive({
//   image: i_norm,
//   description: 'PosthurricaneImage',
//   folder:'PRlndst8',
//   scale: 30,
//   region: studyarea
// });

//Test SMA function try to find the appropriate endmembers
var smalndst_1 = smaima(i_norm);
var smalndst_2 = smaima(i_base);

// Export.image.toDrive({
//   image: smalndst_1,
//   description: 'SMA_Pre',
//   scale: 30,
//   region: PRarea
// });

// Export.image.toDrive({
//   image: smalndst_2,
//   description: 'SMA_Post',
//   scale: 30,
//   region: PRarea
// });

var npv_1 = ee.Image(smalndst_1.select(0));
var gv_1  = ee.Image(smalndst_1.select(1));
var nrmlzd_npv_1 =npv_1.divide(npv_1.add(gv_1));


var npv_2 = ee.Image(smalndst_2.select(0));
var gv_2  = ee.Image(smalndst_2.select(1));
var nrmlzd_npv_2 =npv_2.divide(npv_2.add(gv_2));


var DNPV = npv_1.subtract(npv_2);
var nrmlzdDNPV = nrmlzd_npv_1.subtract(nrmlzd_npv_2);

//Apply a forest mask on DNPV
var hansenforest = hansen.select('treecover2000').clip(PRboundary).gte(50);
var DNPVforest = DNPV.updateMask(hansenforest)
var nrmlzdDNPVforest = nrmlzdDNPV.updateMask(hansenforest)

//Mask exteme values (>1 or <0) in normalized DNPVforest
var gt0 = nrmlzdDNPVforest.gte(0)
var gt1 = nrmlzdDNPVforest.lte(1)
// var nrmlzdDNPVforestEE = nrmlzdDNPVforest.updateMask(gt0).updateMask(gt1);
var nrmlzdDNPVforestEE = nrmlzdDNPVforest;
var nrmlzd_gv_2 = gv_2.divide(npv_2.add(gv_2)).mask(hansenforest).updateMask(gt0).updateMask(gt1)

// Export.image.toDrive({
//   image: nrmlzd_gv_2,
//   description: 'nrmlzdGV2016reTI',
//   folder:'PRlayersPRRegionfromwind',
//   scale: 30,
//   region: studyarea
// });
// Export.image.toDrive({
//   image: nrmlzdDNPVforestEE,
//   description: 'nrmlzdDNPVforestEEreTI',
//   folder:'PRlayersPRRegionfromwind',
//   scale: 30,
//   region: studyarea
// });

Export.image.toAsset({
  image: nrmlzdDNPVforestEE,
  description: 'nrmlzdDNPV_UI',
  scale: 30,
  region: studyarea
});
//Be careful use the intervals, value not in the intervals will be masked
var DNPV_intervals =
  '<RasterSymbolizer>' +
    '<ColorMap  type="intervals" extended="false" >' +
      '<ColorMapEntry color="#4B0082" quantity="0.0" label="-1.0-0.0" />' +        //Indigo
      '<ColorMapEntry color="#4169E1" quantity="0.2" label="0.0001-0.2000" />' +   //royalblue
      '<ColorMapEntry color="#32CD32" quantity="0.4" label="0.2001-0.4000" />' +   //Limegreen
      '<ColorMapEntry color="#FFFF00" quantity="0.6" label="0.4001-0.6000" />' +   //Yellow
      '<ColorMapEntry color="#FFA500" quantity="0.8" label="0.6001-0.8000" />' +   //Orange
      '<ColorMapEntry color="#FF0000" quantity="1.0" label="0.8000-1.0000" />' +   //Red
    '</ColorMap>' +
  '</RasterSymbolizer>';
  
  var DNPVinte =
  '<RasterSymbolizer>' +
    '<ColorMap  type="intervals" extended="false" >' +
      '<ColorMapEntry color="#00FFFF" quantity="-200.0" label="-10.0-0.0" />' +        //black
      '<ColorMapEntry color="#000000" quantity="0.0" label="-1.0-0.0" />' +        //black
      '<ColorMapEntry color="#000000" quantity="0.2" label="0.0001-0.2000" />' +   //black
      '<ColorMapEntry color="#000000" quantity="0.4" label="0.2001-0.4000" />' +   //black
      '<ColorMapEntry color="#000000" quantity="0.6" label="0.4001-0.6000" />' +   //black
      '<ColorMapEntry color="#000000" quantity="1.0" label="0.6001-1.0000" />' +   //black
      '<ColorMapEntry color="#00FFFF" quantity="2.0" label="1.0001-6.5000" />' +   //Red
      '<ColorMapEntry color="#FF0000" quantity="40000.0" label="1.0001-6.5000" />' +   //Red
    '</ColorMap>' +
  '</RasterSymbolizer>';

//Not sure if they use the same standardized symbols
var LULC =
  '<RasterSymbolizer>' +
    ' <ColorMap  type="intervals" extended="false" >' +
    '<ColorMapEntry color="#aec3d4" quantity="0" label="Water"/>' +
    '<ColorMapEntry color="#152106" quantity="1" label="Evergreen Needleleaf Forest"/>' +
    '<ColorMapEntry color="#225129" quantity="2" label="Evergreen Broadleaf Forest"/>' +
    '<ColorMapEntry color="#369b47" quantity="3" label="Deciduous Needleleaf Forest"/>' +
    '<ColorMapEntry color="#30eb5b" quantity="4" label="Deciduous Broadleaf Forest"/>' +
    '<ColorMapEntry color="#387242" quantity="5" label="Mixed Deciduous Forest"/>' +
    '<ColorMapEntry color="#6a2325" quantity="6" label="Closed Shrubland"/>' +
    '<ColorMapEntry color="#c3aa69" quantity="7" label="Open Shrubland"/>' +
    '<ColorMapEntry color="#b76031" quantity="8" label="Woody Savanna"/>' +
    '<ColorMapEntry color="#d9903d" quantity="9" label="Savanna"/>' +
    '<ColorMapEntry color="#91af40" quantity="10" label="Grassland"/>' +
    '<ColorMapEntry color="#111149" quantity="11" label="Permanent Wetland"/>' +
    '<ColorMapEntry color="#cdb33b" quantity="12" label="Cropland"/>' +
    '<ColorMapEntry color="#cc0013" quantity="13" label="Urban"/>' +
    '<ColorMapEntry color="#33280d" quantity="14" label="Crop, Natural Veg. Mosaic"/>' +
    '<ColorMapEntry color="#d7cdcc" quantity="15" label="Permanent Snow, Ice"/>' +
    '<ColorMapEntry color="#f7e084" quantity="16" label="Barren, Desert"/>' +
    '<ColorMapEntry color="#6f6f6f" quantity="17" label="Tundra"/>' +
    '</ColorMap>' +
'</RasterSymbolizer>';

// var DNPVIntensity_intervals =
//   '<RasterSymbolizer>' +
//     '<ColorMap  type="intervals" extended="false" >' +
//       '<ColorMapEntry color="#FFFACD" quantity="-1.0" label="-2.0-0.0" />' +        //lemmonchiffon
//       '<ColorMapEntry color="#FFFACD" quantity="0.0" label="-1.500-0.0000" />' +   //lemonchiffon
//       '<ColorMapEntry color="#F0E68C" quantity="0.2" label="0.0001-0.2000" />' +   //khaki
//       '<ColorMapEntry color="#FFA500" quantity="0.4" label="0.2001-0.4000" />' +   //Orange
//       '<ColorMapEntry color="#FF0000" quantity="0.6" label="0.4001-0.6000" />' +   //Red
//       '<ColorMapEntry color="#A52A2A" quantity="0.8" label="0.6001-0.8000" />' +   //Brown
//       '<ColorMapEntry color="#800080" quantity="1.0" label="0.8001-1.0000" />' +   //Purple
      
//     '</ColorMap>' +
//   '</RasterSymbolizer>';


var DNPVIntensity_intervals =
  '<RasterSymbolizer>' +
    '<ColorMap  type="intervals" extended="false" >' +
      // '<ColorMapEntry color="#FFFACD" quantity="-1.0" label="-2.0-0.0" />' +        //lemmonchiffon
      '<ColorMapEntry color="#442288" quantity="-10.0" label="-0.4000--2.0000" />' + //Spanish Violet
      '<ColorMapEntry color="#6CA2EA" quantity="-0.2" label="-0.2000--0.4000" />' + //Little boy blue
      // '<ColorMapEntry color="#FFFACD" quantity="0.0" label="-1.500-0.0000" />' +   //lemonchiffon
      '<ColorMapEntry color="#B5D33D" quantity="0.0" label="-0.200-0.0000" />' +   //Android Green
      '<ColorMapEntry color="#FFFFE0" quantity="0.2" label="0.0001-0.2000" />' +   //lightyellow
      '<ColorMapEntry color="#FFBB4E" quantity="0.4" label="0.2001-0.4000" />' +   //Pastel Orange
      '<ColorMapEntry color="#FF0000" quantity="0.6" label="0.4001-0.6000" />' +   //Red
      '<ColorMapEntry color="#794B26" quantity="0.8" label="0.6001-0.8000" />' +   //Russet
      '<ColorMapEntry color="#000000" quantity="1.0" label="0.8001-1.0000" />' +   //Black
      
    '</ColorMap>' +
  '</RasterSymbolizer>';  
  
// Map.addLayer(PRLULC.sldStyle(LULC),{}, 'PRLULC')
// Map.addLayer(nrmlzdDNPVforestEE.sldStyle(DNPVinte),{}, 'DNPV larger than 1 in red');
var grey = ee.Image.constant(0).clip(PRboundary)
Map.addLayer(grey, {min:0, max:0, palette:["D3D3D3"]},"grey")
Map.addLayer(nrmlzdDNPVforestEE, {}, 'normalized_DNPV');
Map.addLayer(nrmlzd_npv_2.updateMask(hansenforest).sldStyle(DNPVIntensity_intervals), {}, 'pre_DNPV_clr')
Map.addLayer(nrmlzd_npv_2.updateMask(hansenforest), {}, 'pre_DNPV')
Map.addLayer(nrmlzdDNPVforestEE.sldStyle(DNPVIntensity_intervals), {}, 'normalized_DNPV_clr');


// Export.image.toDrive({
//   image: nrmlzdDNPVforestEE.sldStyle(DNPVIntensity_intervals),
//   description: 'nrmlzdDNPVforestEEinColor',
//   scale: 30,
//   region: studyarea
// });
// // Check the number of DNPV pixels larger than 1
// var output = nrmlzdDNPVforestEE.gt(1)
// var gt1 = nrmlzdDNPVforestEE.updateMask(output);



var point_one = ee.Geometry.Point(-65.763981, 18.311474);
// Map.addLayer(point_one,{},"point 1");

var point_two = ee.Geometry.Point(-66.12739, 18.05068);
// Map.addLayer(point_two,{},"point 2");

var point_three = ee.Geometry.Point(-66.48509, 18.17923);
// Map.addLayer(point_three,{},"point 3");

var point_four = ee.Geometry.Point(-66.69185, 18.32768);
// Map.addLayer(point_four,{},"point 4");

var point_five = ee.Geometry.Point(-65.88425, 18.338372);
// Map.addLayer(point_five,{},"point 5");

var point_six = ee.Geometry.Point(-65.82491, 18.30859)
// Map.addLayer(point_six,{},"point 6")
var elevation = SRTM.clip(PRboundary);
// Map.addLayer(elevation, {},"elevation");

// Map.addLayer(aspect.clip(PRboundary),{},"aspect")
