The purpose of this project is to process a scalar potential map from Comsol to extract a coil winding.

All code is run in python 3, and is probably still compatable with python 2.

Is compatible with Anaconda except for the -s simplification modole of get_field_B0.py requires a non-anaconda package.

This is live code that is still under developement, and comments are sparse in some areas of the code.  Functions are attempted to be named after the operation they perform, but their exact implementation and usage may not be initially obvious.  some legacy code is present that is not used presently is hasn't beem removed yet.

PlottingTools.py - a series of plotting functions used for decorating plots in ways that are common to the application.

get_field_B0.py - main body of the code that loads and interprets the scalar potential file to create a coil model, and then calculates the field  and coil metrics. Most options for how to run the code are set through command line interface options

patch.py - uses the Biot-Savart Law to calculate the fields from a coil and includes from plotting functions for visuallizing a coil.

PipesFitting.py - given a CoilSet class object from patch.py, re-routes the coils around the cylindrical pipes defined by the user. Note: wires can only come in on a principle axis, the pipes must know the direction before running, and no points in the coil to be rerouted can be within the pipe.  These simplifications were made to reduce developement time.

get_field_InputFiles - folder for input files to get_field_B0.  One file is included.

get_field_OutputFiles - folder for the outputs from get_field_B0.py

B0-WireswithShielding.mph - example of using a processed contours from get_field_B0 to simiulate the effect of a surrounding mu-metal shield on the field produced by the coil.  The WiresLoad application builder method is used to load the wire paths from the csv files produced get_field_b0.  Due to GitHub drive space contraints all but one wire path were removed from the model, and the csv files need to be regenerated before it can be run.s


Example code execution:

To read in a data file of scalar potential maps, chosing the contour locations from x.txt. -gtc sets a number of useful plots to be made
python get_field_B0 -f get_field_InputFiles/data-2.0.txt -x -gtc

Performs the optimization of the choice of the scalar potential levels starting with the values given in x.txt
python get_field_B0 -f get_field_InputFiles/data-2.0.txt -x -o

Performs the optimization of the choice of the scalar potential levels starting with an equal spacing of 0.17amps.
python get_field_B0 -f get_field_InputFiles/data-2.0.txt -i 0.17 -o
