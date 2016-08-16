#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
#  Copyright Kitware Inc. and Epidemico Inc.
#
#  Licensed under the Apache License, Version 2.0 ( the "License" );
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
##############################################################################
import gdal
import os
import ogr
import numpy as np
import itertools
from gaia.gaia_process import GaiaProcess
from gaia.geo.geo_inputs import VectorFileIO
from skimage.graph import route_through_array
from math import sqrt, ceil
try:
    import osr
except ImportError:
    from osgeo import osr
import gaia.formats as formats


class LeastCostProcess(GaiaProcess):
    """
    Least cost path analysis.
    """

    default_output = formats.JSON

    def __init__(self, **kwargs):
        super(LeastCostProcess, self).__init__(**kwargs)

        if not self.output:
            self.output = VectorFileIO(
                name='result',
                uri=self.get_outpath())
        self.validate()

        if self.inputs:
            self.start_point = self.inputs[0]['start']
            self.end_point = self.inputs[0]['end']
            self.raster_layer = self.inputs[0]['uri']

    def pixel_offset2coord(self, rasterfn, xOffset, yOffset):
        raster = gdal.Open(rasterfn)
        geotransform = raster.GetGeoTransform()
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]
        coordX = originX+pixelWidth*xOffset
        coordY = originY+pixelHeight*yOffset
        return coordX, coordY

    def array2shp(self, array, outSHPfn, rasterfn, pixelValue):
        raster = gdal.Open(rasterfn)
        geotransform = raster.GetGeoTransform()
        pixelWidth = geotransform[1]
        maxDistance = ceil(sqrt(2*pixelWidth*pixelWidth))

        count = 0
        roadList = np.where(array == pixelValue)
        pointDict = {}
        for indexY in roadList[0]:
            indexX = roadList[1][count]
            Xcoord, Ycoord = self.pixel_offset2coord(rasterfn, indexX, indexY)
            pointDict[count] = (Xcoord, Ycoord)
            count += 1

        multiline = ogr.Geometry(ogr.wkbMultiLineString)
        for i in itertools.combinations(pointDict.values(), 2):
            point1 = ogr.Geometry(ogr.wkbPoint)
            point1.AddPoint(i[0][0], i[0][1])
            point2 = ogr.Geometry(ogr.wkbPoint)
            point2.AddPoint(i[1][0], i[1][1])

            distance = point1.Distance(point2)

            # calculate the distance between two points
            if distance < maxDistance:
                line = ogr.Geometry(ogr.wkbLineString)
                line.AddPoint(i[0][0], i[0][1])
                line.AddPoint(i[1][0], i[1][1])
                multiline.AddGeometry(line)

        shpDriver = ogr.GetDriverByName("GeoJSON")
        if os.path.exists(outSHPfn):
            shpDriver.DeleteDataSource(outSHPfn)
        else:
            self.output.create_output_dir(outSHPfn)

        outDataSource = shpDriver.CreateDataSource(outSHPfn)
        outLayer = outDataSource.CreateLayer(outSHPfn,
                                             geom_type=ogr.wkbMultiLineString)

        featureDefn = outLayer.GetLayerDefn()
        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(multiline)
        outLayer.CreateFeature(outFeature)

    def raster_to_array(self, raster):
        raster = gdal.Open(raster)
        band = raster.GetRasterBand(1)
        array = band.ReadAsArray()
        return array

    def coord2pixeloffset(self, rasterfn, x, y):
        raster = gdal.Open(rasterfn)
        geotransform = raster.GetGeoTransform()
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]
        xOffset = int((x - originX)/pixelWidth)
        yOffset = int((y - originY)/pixelHeight)
        return xOffset, yOffset

    def create_path(self, raster, costSurfaceArray, start, end):

        # coordinates to array index
        startCoordX = start[0]
        startCoordY = start[1]
        startIndexX, startIndexY = self.coord2pixeloffset(
            raster, startCoordX, startCoordY)

        stopCoordX = end[0]
        stopCoordY = end[1]
        stopIndexX, stopIndexY = self.coord2pixeloffset(
            raster, stopCoordX, stopCoordY)

        # create path
        indices, weight = route_through_array(
            costSurfaceArray,
            (startIndexY, startIndexX),
            (stopIndexY, stopIndexX),
            geometric=True,
            fully_connected=True)
        indices = np.array(indices).T
        path = np.zeros_like(costSurfaceArray)
        path[indices[0], indices[1]] = 1
        return path

    def array_to_raster(self, newRasterfn, rasterfn, array):
        raster = gdal.Open(rasterfn)
        geotransform = raster.GetGeoTransform()
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]
        cols = array.shape[1]
        rows = array.shape[0]

        driver = gdal.GetDriverByName('GTiff')
        outRaster = driver.Create(newRasterfn, cols, rows, gdal.GDT_Byte)
        outRaster.SetGeoTransform(
            (originX, pixelWidth, 0, originY, 0, pixelHeight))
        outband = outRaster.GetRasterBand(1)
        outband.WriteArray(array)
        outRasterSRS = osr.SpatialReference()
        outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()

    def calculate_path(self, raster, start, end):
        costSurfaceArray = self.raster_to_array(raster)
        pathArray = self.create_path(raster, costSurfaceArray, start, end)
        self.array2shp(pathArray, self.output.uri, raster,  1)

    def compute(self):
        self.calculate_path(
            self.raster_layer, self.start_point, self.end_point)


PLUGIN_CLASS_EXPORTS = [
    LeastCostProcess,
]
