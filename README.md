# AMTT2.0
Automated Markland's Test Tool, v 2.0
The code refers to an automatic procedure based on the GIS environment working principles and developed it in Matlab language. 
The Matlab code performs Markland's tests for planar and wedge sliding, direct and flexural toppling.
With respect to the original version of AMTT, this new one allows the user to choose among two different modes for performing the kinematic tests, depending on the input data regarding the discontinuities: the Deterministic mode and the Raw Data mode.
The Deterministic mode requires the orientation of the main sets, the Raw Data mode is capable of considering directly the raw data from a geomechanical survey.
For both modes, discontinuities orientation and friction angles, along with slope and aspect data representing the rockface orientation of the considered outcrop, are the input data. 
The slope and aspect data must be in TIFF or GeoTIFF format. 

The outputs are a series of GeoTIFF raster files describing the result for each kinematism separately and globally, which can be imported directly into GIS, with the same extent and georeferencing of the input data.
The global results can be also used to map source areas for 3D rockfall numerical simulations. 

How to cite the original code and the new version:

Taboni, B.; Tagliaferri, I.D.; Umili, G. A Tool for Performing Automatic Kinematic Analysis on Rock Outcrops. Geosciences 2022, 12, 435. https://doi.org/10.3390/geosciences12120435 
