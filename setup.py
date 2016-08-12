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
from setuptools import setup, find_packages

setup(
  name="gaia-leastcostpath-plugin",
  version="0.0",
  description="""Gaia plugin""",
  author="Essam",
  install_requires=["gaia>=0.0.0"],
  packages=find_packages(),
  include_package_data=True,
  entry_points={
    'gaia.plugins': [
            "gaia_leastcostpath.processes = gaia_leastcostpath.processes",
        ]
  }
)
