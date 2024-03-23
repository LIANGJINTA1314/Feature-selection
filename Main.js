/*
"roi represents the selected study area. 
'construction', 'ocean', 'water', 'wood', 'field', 'bare', 'tidalflat', 'wetland' are selected sample points."
*/

//Inversion of LST
    var landsat8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")

    var LST = landsat8.filterDate('2021-01-01', '2021-12-31')
                            .filterBounds(roi)
                            .filter(ee.Filter.lte('CLOUD_COVER',10))
                            .select('SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','ST_B10','QA_PIXEL','QA_RADSAT')
                            .sort('system:index')
                            .map(CalLST)
                            .median()
                            .clip(roi)

    function CalLST(image) {
      var b1 = image.select('ST_B10');
      var Albedo = image.select('ST_B10');
      return Albedo.multiply(0.00341802).add(149.0).subtract(273.15).rename(['LST']).copyProperties(image, image.propertyNames());
    }

    Map.addLayer(LST, {'min':2,'max':49,'palette':["eff3ff","c6dbef","9ecae1","6baed6","4292c6","2171b5","084594",
    "fff5f0","fee0d2","fcbba1","fc9272","fb6a4a","ef3b2c","cb181d","99000d"]}, 'LST',false)


//Calculation of slope and aspect
    var srtm = ee.Image('USGS/SRTMGL1_003');
    var elevation = srtm.select('elevation').clip(roi);

    var slope = ee.Terrain.slope(elevation).clip(roi);

    var aspect = ee.Terrain.aspect(srtm);
    // Stretch the slope aspect values to be between -1 and 1
    //-1 represents due west, and 1 represents due east direction
    var asImage = aspect.divide(180).multiply(Math.PI).sin().clip(roi);
    Map.addLayer(slope, {min: 0, max:30}, 'slope',false);
    Map.addLayer(asImage, {min:-1,max:1}, 'asImage',false)


//Import components after Principal Component Analysis (PCA)
    var pcaVisParams = {bands: ['pc1', 'pc2', 'pc3'], min: -2, max: 2};
    Map.addLayer(pcaImported, pcaVisParams, 'Principal Components',false);


//Calculation of various indices
    //SAVI
    function S2_savi(image) {
      // Add Soil Adjust Vegetation Index (SAVI)
        // using L = 0.5;
        var savi = image.expression('(NIR - RED) * (1 + 0.5)/(NIR + RED + 0.5)', {
        'NIR': image.select('B8'),
        'RED': image.select('B4')
        }).float();
        return image.addBands(savi.rename('SAVI'));
    }

    //EVI
    function S2_evi(image){
         var evi = image.expression('2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
        'NIR' : image.select('B8'),
        'RED' : image.select('B4'),
        'BLUE': image.select('B2')
      }).float();
      return image.addBands(evi.rename('EVI'));
    }

    // RVI
    function S2_rvi(image){
        var rvi = image.expression('NIR / Red', {
        'NIR': image.select('B8'),
        'Red': image.select('B4')
      }).float();
      return image.addBands(rvi.rename('RVI'));
     }


//Add features TIR1 and TIR2 that are not present in Sentinel.
    //Landsat 8 image quality screening mask
    function maskL8sr(image) {
      // The third and fifth bands respectively represent cloud shadow and cloud
      var cloudShadowBitMask = 1 << 3;
      var cloudsBitMask = 1 << 5;
 
      // Retrieve the pixel QA band
      var qa = image.select('pixel_qa');
 
      // "Set both values to 0 under specific conditions."
      var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
          .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
 
      // Update the cloud mask bands, then scale by reflectance, select band attributes, and finally assign to the image
      return image.updateMask(mask).divide(10000)
          .select("B[0-9]*")
          .copyProperties(image, ["system:time_start"]);
    }

    var L8_SR = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
        .filterDate('2020-01-01', '2020-12-31')
        .map(maskL8sr)
        .median()
    var BL10 = L8_SR.select('B10')
    var BL11 = L8_SR.select('B11')


//Sentinel-2 image quality screening mask
    function maskS2clouds(image) {
      var qa = image.select('QA60');

      // Bits 10 and 11 are clouds and cirrus, respectively.
      var cloudBitMask = 1 << 10;
      var cirrusBitMask = 1 << 11;

      // Both flags should be set to zero, indicating clear conditions.
      var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
          .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

      return image.updateMask(mask).divide(10000);
    }


//Create annual composite image
    // Load Sentinel-2 SR reflectance data.
    var collection = ee.ImageCollection('COPERNICUS/S2_SR')
                  .filterBounds(roi)
                  .filterDate('2021-01-01', '2021-12-31')
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10))
                  .map(maskS2clouds)
                  .map(S2_savi)
                  .map(S2_evi)
                  .map(S2_rvi)
                  
    //Assign median to each pixel, clip the image to the study area
    var image1 = collection.median().clip(roi);
    
    // Display the image with 12 11 4
    var rgbVis = {
      min: 0.0,
      max: 0.3,
      bands: ['B12', 'B11', 'B4'],
    };
    Map.addLayer(image1, rgbVis, '12 11 4',false);

    // Display the image with 4 3 2
    Map.addLayer(image1, {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3},'4 3 2',false);


//Calculate the remaining indices
    var mndwi = image1.normalizedDifference(['B3', 'B11']).rename('MNDWI');
    var ndbi = image1.normalizedDifference(['B11', 'B8']).rename('NDBI');
    var ndvi = image1.normalizedDifference(['B8', 'B4']).rename('NDVI');


//Extract the first three components from the image after Principal Component Analysis (PCA)
    var pc1 = pcaImported.select("pc1")
    var pc2 = pcaImported.select("pc2")
    var pc3 = pcaImported.select("pc3")


//Compute glcm texture feature
    var gray = image1.expression(
          '(0.3 * NIR) + (0.59 * R) + (0.11 * G)', {
          'NIR': image1.select('B8'),
          'R': image1.select('B4'),
          'G': image1.select('B3'),
    }).rename('gray');
 
    var glcm = gray.unitScale(0,0.30).multiply(100).toInt().glcmTexture({size: 1,kernel:null});
    print('glcm',glcm)

    //Extract certain features and rename them
    var SAVG = glcm.select('gray_savg').rename('SAVG');
    var CORR = glcm.select('gray_corr').rename('CORR');
    var VAR = glcm.select('gray_var').rename('VAR');
    var DENT = glcm.select('gray_dent').rename('DENT');
    var ASM = glcm.select('gray_asm').rename('ASM');
    var CONT = glcm.select('gray_contrast').rename('CONT');
    var IDM = glcm.select('gray_idm').rename('IDM');
    var SENT = glcm.select('gray_sent').rename('SENT');
    var ENT = glcm.select('gray_ent').rename('ENT');


//Add features into image
    image1=image1.addBands(ndvi).addBands(ndbi).addBands(mndwi)
    .addBands(SAVG).addBands(CORR).addBands(VAR).addBands(DENT).addBands(CONT).addBands(IDM).addBands(ASM)
    .addBands(SENT).addBands(ENT).addBands(elevation.rename('DEM'))
    .addBands(pc1).addBands(pc2).addBands(pc3)
    .addBands(slope.rename('Slope')).addBands(asImage.rename('Aspect'))
    .addBands(LST).addBands(BL10.rename('BL10')).addBands(BL11.rename('BL11'))
    image1=image1.clip(roi)
    
    //The current number of features being inputted here
    var bands1 = ['B2', 'B3', 'B4', 'B8','B11','B12','EVI','SAVI','MNDWI','NDBI','NDVI','SAVG',
    'CORR','VAR','DENT','pc1','pc2','pc3','Slope','Aspect','LST','RVI','DEM','ASM','CONT','IDM',
    'SENT','ENT','BL10','BL11'];


//Integrate sample points
    var classNames1 = construction.merge(oceam).merge(water).merge(wood).merge(field).merge(bare).merge(tidalflat).merge(wetland);
    var training1 = image1.select(bands1).sampleRegions({
      collection: classNames1,
      properties: ['landcover'],
      scale: 10,
      tileScale:8
    });


//random uniforms to the training dataset.
    var withRandom1 = training1.randomColumn('random');


//Retain some data for testing to prevent overfitting of the model
    var split = 0.7; 
    var trainingPartition1 = withRandom1.filter(ee.Filter.lt('random', split));
    var testingPartition1 = withRandom1.filter(ee.Filter.gte('random', split));


//Select attributes for classification
    var classProperty = 'landcover';



//Random forest parameter settings
    var numberOfTrees = 200; //200-600
    var variablesPerSplit = null; 
    var minLeafPopulation = 1; 
    var bagFraction = 0.5; 
    var maxNodes = null; 
    var seed = 0;
    
    var classifier1 = ee.Classifier.smileRandomForest({
      numberOfTrees: numberOfTrees,
      variablesPerSplit: variablesPerSplit,
      minLeafPopulation: minLeafPopulation,
      bagFraction: bagFraction,
      maxNodes: maxNodes,
      seed: seed
    }).train({
      features: trainingPartition1,
      classProperty: 'landcover',
      inputProperties: bands1,
    });

    var classified1 = image1.select(bands1).classify(classifier1);


//Display the classification map
    Map.centerObject(classNames1, 11);
   
   //Color settings
    var imageVisParam1 ={
      min:1,
      max:8,
      opacity:1,
      palette:["#FF0000","#0000FF","#19c9f9","#138d0f","#54ca0b","#b19c07","#fbff04","#248684"]
    } 
    
    Map.addLayer(classified1,imageVisParam1,'classified1');

//Feature importance
    var dict1 = classifier1.explain();
    print('Explain1:',dict1);
     
    var variable_importance1 = ee.Feature(null, ee.Dictionary(dict1).get('importance'));
    var chart1 =
    ui.Chart.feature.byProperty(variable_importance1)
    .setChartType('ColumnChart')
    .setOptions({
    title: 'Random Forest Variable Importance',
    legend: {position: 'none'},
    hAxis: {title: 'Bands'},
    vAxis: {title: 'Importance'}
    });
     
    print(chart1);


//Confusion matrix and accuracy validation
    var test1 = testingPartition1.classify(classifier1);
    
    var confusionMatrix1 = test1.errorMatrix('landcover', 'classification');
    print('confusionMatrix1',confusionMatrix1);
    print('overall accuracy1', confusionMatrix1.accuracy());
    print('kappa accuracy1', confusionMatrix1.kappa());


//Input classification results to Drive
    var img1 = classified1.visualize(imageVisParam1)
    
    Export.image.toDrive({
      image:img1,
      description:'jiangsu',
      region:roi,
      scale:10,
      maxPixels:1e13
    })



