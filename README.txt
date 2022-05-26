The purpose of this project is to process a scalar potential map from Comsol to extract a coil winding.

This is live code that is still under developement, and comments are sparse in some areas of the code.  Functions are attempted to be named after the operation they perform, but their exact implementation and usage may not be initially obvious.  some legacy code is present that is not used presently is still present in the code.

PlottingTools.py - a series of plotting functions used for decorating plots in ways that are common to the application.

get_field_B0.py - main body of the code that loads and interprets the scalar potential file to create a coil model, and then calculates the field  and coil metrics.

patch.py - uses the Biot-Savart Law to calculate the fields from a coil and includes from plotting functions for visuallizing a coil.

PipesFitting.py - given a CoilSet class object from patch.py, re-routes the coils around the cylindrical pipes defined by the user. Note: wires can only come in on a principle axis, the pipes must know the direction before running, and no points in the coil to be rerouted can be within the pipe.  These simplifications were made to reduce developement time.

get_field_InputFiles - folder for input files to get_field_B0.  One file is included.

get_field_OutputFiles - folder for the outputs from get_field_B0.py