## gaia-least-cost-path-plugin

This is a plugin for Gaia (https://github.com/OpenDataAnalytics/gaia) that
computes the 'least cost path' between two points based on a grid of cost values.
Least cost path analysis calculates the most cost-effective route between a source and destination.
Cost can be a function of elevation, time, or any other criteria that is represented as values on a raster grid, where a higher value indicates a higher cost.
As part of the analysis, the 8 neighbors of a grid cell are evaluated and the path moves to the cell with the smallest value.
This evaluation is repeated until the source and destination are connected.  The output is a vector line that connects the source and destination points.


#### Documentation

Documentation for Gaia can be found at http://gaia.readthedocs.org.
Documentation for this plugin can be found at http://gaia-leastcostpath-plugin.readthedocs.org.

#### Installation

  - pip install -e .
  - pip install -r requirements.txt

#### Inputs
  - A raster image
  - Start and end point

#### License

Copyright 2015 Kitware Inc.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0


Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
