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
from gaia.parser import deserialize
from gaia.core import config

testfile_path = os.path.join(os.path.dirname(
    os.path.realpath(__file__)), '../data')


class TestLeastCostViaParser(unittest.TestCase):
    """Tests for the Gaia Least Cost plugin via Parser"""

    def test_process_least_cost_path_2vectors(self):
        """Test Least Cost Process"""
        with open(os.path.join(testfile_path,
                               'least_cost_path_2vectors.json')) as inf:
            body_text = inf.read().replace('{basepath}', testfile_path)
        process = json.loads(body_text, object_hook=deserialize)
        try:
            process.compute()
            output = json.loads(process.output.read(format=formats.JSON))
            with open(os.path.join(
                    testfile_path,
                    'least_cost_path_process_results.json')) as gj:
                expected_json = json.load(gj)
            self.assertIn('features', output)

            self.assertEquals(len(expected_json['features']),
                              len(output['features']))
            self.assertIsNotNone(process.id)
            self.assertIn(process.id, process.output.uri)
        finally:
            if process:
                process.purge()


def test_process_least_cost_path_1vector_2points(self):
    """Test Least Cost Process"""
    with open(os.path.join(testfile_path,
                           'least_cost_path_1vector.json')) as inf:
        body_text = inf.read().replace('{basepath}', testfile_path)
    process = json.loads(body_text, object_hook=deserialize)
    try:
        process.compute()
        output = json.loads(process.output.read(format=formats.JSON))
        with open(os.path.join(
                testfile_path,
                'least_cost_path_process_results.json')) as gj:
            expected_json = json.load(gj)
        self.assertIn('features', output)

        self.assertEquals(len(expected_json['features']),
                          len(output['features']))
        self.assertIsNotNone(process.id)
        self.assertIn(process.id, process.output.uri)
    finally:
        if process:
            process.purge()