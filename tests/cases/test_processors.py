#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
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
###############################################################################
import json
import os
import unittest
from gaia import formats
from gaia.geo import RasterFileIO
from gaia.geo.geo_inputs import FeatureIO
from gaia_leastcostpath.processes import LeastCostProcess

testfile_path = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), '../data')


class TestLeastCostProcessors(unittest.TestCase):
    def test_process_least_cost_path(self):
        """
        Test LeastCostProcess for raster inputs
        """
        raster = RasterFileIO(
            uri=os.path.join(testfile_path, 'globalprecip.tif'))
        vector1 = FeatureIO(
            features=[{"geometry":
                      {"type": "Point",
                       "coordinates": [-71.590526, 42.659566]
                       },
                       "properties": {}
                       }]
        )
        vector2 = FeatureIO(
            features=[{"geometry":
                      {"type": "Point",
                       "coordinates": [-122.817426, 46.561380]},
                       "properties": {}
                       }]
        )

        process = LeastCostProcess(inputs=[raster, vector1, vector2])
        try:
            process.compute()
            with open(os.path.join(
                    testfile_path,
                    'least_cost_path_process_results.json')) as exp:
                expected_json = json.load(exp)
            actual_json = json.loads(process.output.read(format=formats.JSON))
            self.assertEquals(len(expected_json['features']),
                              len(actual_json['features']))
        finally:
            if process:
                process.purge()
