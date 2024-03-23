/*
This code performs Principal Component Analysis (PCA) on remote sensing imagery 
to generate 'pcaImported' required in the 'Main'.
*/


// Sentinel-2 image quality screening mask
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


// Map the function over one year of data.
// Load Sentinel-2 SR reflectance data.
    var dataset = ee.ImageCollection('COPERNICUS/S2_SR')
                      .filterBounds(roi)  // your study area
                      .filterDate('2021-01-01', '2021-12-31')
                      // Pre-filter to get less cloudy granules.
                      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 10))
                      .map(maskS2clouds);
    
    var s2 = dataset.median().clip(roi);


// A function to calculate principal components
function PCA(maskedImage) {
  var image = maskedImage.unmask();
  var scale = 20;
  var region = roi;
  var bandNames = image.bandNames();

  // Perform mean and median processing on the data to achieve faster covariance statistical calculations and standard deviation stretching of principal components
  var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9,
    bestEffort: true,
    tileScale: 16
  });
  var means = ee.Image.constant(meanDict.values(bandNames));
  var centered = image.subtract(means);

  // This helper function returns a list of new band names
  var getNewBandNames = function(prefix) {
    var seq = ee.List.sequence(1, bandNames.length());
    return seq.map(function(b) {
      return ee.String(prefix).cat(ee.Number(b).int());
    });
  };

  // This function accepts an image with mean centering, a scale factor, and a region for analysis 
  // It returns the principal components (PCs) within the specified region as a new image
  var getPrincipalComponents = function(centered, scale, region) {
    // Flatten the bands of the image into a one-dimensional array for each pixel
    var arrays = centered.toArray();
    
    // Calculate the covariance of the bands within the region
    var covar = arrays.reduceRegion({
      reducer: ee.Reducer.centeredCovariance(),
      geometry: region,
      scale: scale,
      maxPixels: 1e9,
      bestEffort: true,
      tileScale: 16
    });

    // Obtain the array covariance results and convert them into a single array
    // This represents the covariance between bands within the region
    var covarArray = ee.Array(covar.get('array'));
    
    // Perform feature analysis, separating numerical and vector components
    var eigens = covarArray.eigen();
 
    // Eigenvalue vector
    var eigenValues = eigens.slice(1, 0, 1);
    
    // Calculate the percentage difference for each component
    var eigenValuesList = eigenValues.toList().flatten();
    var total = eigenValuesList.reduce(ee.Reducer.sum());
    var percentageVariance = eigenValuesList.map(function(item) {
      return (ee.Number(item).divide(total)).multiply(100).format('%.2f');
    });

    // This will allow us to determine how many components capture most of the variance in the input
    print('Percentage Variance of Each Component', percentageVariance);

    // This is a PxP matrix where the rows contain the eigenvectors。
    var eigenVectors = eigens.slice(1, 1);

    // Convert the array image into a two-dimensional array for matrix calculations
    var arrayImage = arrays.toArray(1);

    // Left-multiply the image array by the eigenvector matrix
    var principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage);

    // Square root of eigenvalues turned into P-band image
    var sdImage = ee.Image(eigenValues.sqrt())
      .arrayProject([0]).arrayFlatten([getNewBandNames('sd')]);

    // Convert principal components into P-band image, normalized by standard deviation (SD)
    return principalComponents
      // Drop an unnecessary dimension, [[]]. -> [].
      .arrayProject([0])
      // Convert an array image of one band into a multi-band image. [] -> Image
      .arrayFlatten([getNewBandNames('pc')])
      // Normalize the principal components according to their standard deviations
      .divide(sdImage);
  };

  var pcImage = getPrincipalComponents(centered, scale, region);
  return pcImage.mask(maskedImage.mask());
}

// Obtain pcaImported
var pcaImported = PCA(s2).select(['pc1', 'pc2', 'pc3']);
print(pcaImported);


