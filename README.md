# GEE_rsei
This is a code for calculating the remote sensing ecological index, how should I modify it to download the image of four components such as NDVI
var roi = ee.FeatureCollection(table)
var star_date = '2018-01-01'//定义起始时间
var end_date = '2018-12-30'//定义终止时间
var L8_SR = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")//加载L8_SR影像
var img_SR = L8_SR.filterBounds(roi)
           .filterDate(star_date, end_date)
           .filterMetadata('CLOUD_COVER', 'less_than',50)
          .mean()//按照条件进行筛选
           
//归一化函数
var img_normalize = function(img){
      var minMax = img.reduceRegion({
            reducer:ee.Reducer.minMax(),
            geometry: roi,
            scale: 100,
            maxPixels: 10e13,
            tileScale: 16
        })
      var year = img.get('year')
      var normalize  = ee.ImageCollection.fromImages(
            img.bandNames().map(function(name){
                  name = ee.String(name);
                  var band = img.select(name);
                  return band.unitScale(ee.Number(minMax.get(name.cat('_min'))), ee.Number(minMax.get(name.cat('_max'))));
                    
              })
        ).toBands().rename(img.bandNames());
        return normalize;
} 

// Load the image of landsat 8 data

var vizParams = {bands: ['B5', 'B4', 'B3'], min: 50, max: 1500, gamma: [0.95,1.1,1]};

var images = img_SR.clip(roi);

var thermal = images.select('B10').multiply(0.1)

var Emissivity = function(image, nir, red, min, max, de) {
  min = min || 0.2
  max = max || 0.5
  de = de || 0.005
  
  var ndvi_e = "(nir-red)/(nir+red)"
  var ndvi =  ee.Image(0).expression(ndvi_e, {'nir': image.select(nir), 'red': image.select(red)}).rename('ndvi')
  
  var pv = ee.Image(0).expression("((ndvi-min)/(max-min))**2", {ndvi: ndvi, min: min, max: max})
  
  var exp = "ndvi < 0.2 ? 0.979 - (0.046 * b4) : 0.2 <= ndvi <= 0.5 ? (0.987 * pv) + 0.971 * (1 - pv) + de : 0.987 + de"
  
  return ee.Image(0).expression(exp, {ndvi: ndvi, b4: image.select(red).multiply(0.0001), pv: pv, de: de}).rename('EM')
}

var EMM1 = Emissivity(img_SR, 'B5', 'B4')
var EMM = EMM1.clip(roi);



//LST计算
var LST = thermal.expression(
    '(Tb/(1 + (0.0010904 * (Tb / 1.438))*log(Ep)))-273.15', {
      'Tb': thermal.select('B10'),
      'Ep': EMM.select('EM'),
     
});
 LST = img_normalize(LST)
img_SR = img_SR.addBands(LST.rename('LST').toFloat())
var SR_LST= img_SR.select('LST')

// 去云函数 
function removeCloud(image){
  var qa = image.select('BQA')
  var cloudMask = qa.bitwiseAnd(1 << 4).eq(0)
  var cloudShadowMask = qa.bitwiseAnd(1 << 8).eq(0)
  var valid = cloudMask.and(cloudShadowMask)
  return image.updateMask(valid)
}
var L8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA')
var img = L8.filterBounds(roi)
           .filterDate(star_date, end_date)
           .filterMetadata('CLOUD_COVER', 'less_than',50)
           .map(removeCloud).mean()
           
//计算NDVI

var ndvi = img.normalizedDifference(['B5', 'B4']);
  img = img.addBands(ndvi.rename('NDVI'))

//计算WET

//WET
var Wet = img.expression('B*(0.1509) + G*(0.1973) + R*(0.3279) + NIR*(0.3406) + SWIR1*(-0.7112) + SWIR2*(-0.4572)',{
       'B': img.select(['B2']),
       'G': img.select(['B3']),
       'R': img.select(['B4']),
       'NIR': img.select(['B5']),
       'SWIR1': img.select(['B6']),
       'SWIR2': img.select(['B7'])
     })   
  Wet = img_normalize(Wet)
  img = img.addBands(Wet.rename('WET').toFloat())
  
  //计算NDBSI
  
  var ibi = img.expression('(2 * SWIR1 / (SWIR1 + NIR) - (NIR / (NIR + RED) + GREEN / (GREEN + SWIR1))) / (2 * SWIR1 / (SWIR1 + NIR) + (NIR / (NIR + RED) + GREEN / (GREEN + SWIR1)))', {
      'SWIR1': img.select('B6'),
      'NIR': img.select('B5'),
      'RED': img.select('B4'),
      'GREEN': img.select('B3')
    })
   
  var si = img.expression('((SWIR1 + RED) - (NIR + BLUE)) / ((SWIR1 + RED) + (NIR + BLUE))', {
      'SWIR1': img.select('B6'),
      'NIR': img.select('B5'),
      'RED': img.select('B4'),
      'BLUE': img.select('B2')
    }) 
  var ndbsi = (ibi.add(si)).divide(2)
   ndbsi= img_normalize(ndbsi)
  img = img.addBands(ndbsi.rename('NDBSI'))
  var L8_img = img.select(["NDVI","WET","NDBSI"])
  L8_img = SR_LST.addBands(L8_img)

  
//输入基本信息 

var region = roi;
var image =  sentImage.select(bands);
var scale = 100;
var bandNames = image.bandNames();

// 图像波段重命名函数
var getNewBandNames = function(prefix) {
    var seq = ee.List.sequence(1, bandNames.length());
    return seq.map(function(b) {
      return ee.String(prefix).cat(ee.Number(b).int());
    });
  };
//数据平均

var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);


//主成分分析函数
var getPrincipalComponents = function(centered, scale, region) {
    // 图像转为一维数组
    var arrays = centered.toArray();

    // 计算相关系数矩阵
    var covar = arrays.reduceRegion({
      reducer: ee.Reducer.centeredCovariance(),
      geometry: region,
      scale: scale,
      maxPixels: 1e9
    });

    // 获取“数组”协方差结果并转换为数组。
    // 波段与波段之间的协方差
    var covarArray = ee.Array(covar.get('array'));

    // 执行特征分析，并分割值和向量。
    var eigens = covarArray.eigen();

    // 特征值的P向量长度
    var eigenValues = eigens.slice(1, 0, 1);
   

    //计算主成分载荷
    var eigenValuesList = eigenValues.toList().flatten()
    var total = eigenValuesList.reduce(ee.Reducer.sum())
    var percentageVariance = eigenValuesList.map(function(item) {
      return (ee.Number(item).divide(total)).multiply(100).format('%.2f')
    })
 print(eigenValues ,'特征值')
    print("贡献率", percentageVariance)  

    // PxP矩阵，其特征向量为行。
    var eigenVectors = eigens.slice(1, 1);

    // 将图像转换为二维阵列
    var arrayImage = arrays.toArray(1);

    //使用特征向量矩阵左乘图像阵列
    var principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage);

    // 将特征值的平方根转换为P波段图像。
    var sdImage = ee.Image(eigenValues.sqrt())
      .arrayProject([0]).arrayFlatten([getNewBandNames('sd')]);

    //将PC转换为P波段图像，通过SD标准化。
    principalComponents=principalComponents
      // 抛出一个不需要的维度，[[]]->[]。
      .arrayProject([0])
      // 使单波段阵列映像成为多波段映像，[]->image。
      .arrayFlatten([getNewBandNames('pc')])
      // 通过SDs使PC正常化。
      .divide(sdImage);
    return principalComponents
  };
//进行主成分分析，获得分析结果
var pcImage = getPrincipalComponents(centered, scale, region);
//RESI
var img = img.addBands(ee.Image(1).rename('constant'))

var rsei = img.expression('constant - pc1' , {
             constant: img.select('constant'),
             pc1: pcImage.select('pc1')
         })
rsei = img_normalize(rsei)
rsei = pcImage.addBands(rsei.rename('rsei'))



print(rsei,"RSEI")
var visParam = {
    palette: 'FFFFFF, CE7E45, DF923D, F1B555, FCD163, 99B718, 74A901, 66A000, 529400,' +
        '3E8601, 207401, 056201, 004C00, 023B01, 012E01, 011D01, 011301'
 };
//添加图层
 Map.addLayer(rsei.select('rsei').clip(roi), visParam, 'rsei')
//计算RSEI平均值
var rsei_mean = image.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry:roi,
  scale: 30,
  maxPixels: 1e13
});
print(rsei_mean,'平均值')
var rsei_std = image.reduceRegion({
  reducer: ee.Reducer.stdDev(),
  geometry:roi,
  scale: 30,
  maxPixels: 1e13
});
print(rsei_std,'标准差')
// var resimax = resi.select('rsei').reduceRegion({
//   reducer: ee.Reducer.max(),
//   geometry:roi,
//   scale: 30,
//   maxPixels: 1e13
// });
//下载函数
Export.image.toDrive({
   image: rsei.select(["rsei"]),
   description: 'rsei',
   scale: 30,
   region:roi,
   fileFormat: 'GeoTIFF',
 });
