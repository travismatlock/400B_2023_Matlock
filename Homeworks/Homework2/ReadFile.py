# Author: Travis Matlock
# ASTR 400B 1/24/23

# Import modules
import numpy as np
import astropy.units as u

# This function takes the a file (MW_nnn) and reads it out into usable data.
# It takes the name of the file to read out as a string. It returns the time of
# the simulation in Myr, the total number of particles in the simulation, and
# a data array of the type, mass, and components of position and velocity for
# each particle
def Read(filename):
    # Open the file for manipulation
    file = open(filename, 'r')
    # Read the first two lines
    line1 = file.readline()
    line2 = file.readline()
    # Create a list of label and value for each line. Assign variables to each
    # element.
    label_1, value_1 = line1.split()
    label_2, value_2 = line2.split()
    # Convert the time to a float and assign units
    time = float(value_1)*u.Myr
    # Convert the total number of particles to a float
    total = float(value_2)
    file.close()
    # Generate a full data array using the rest of the txt file.
    # dtype=None allows the program to determine the type of each array element
    # regularly (not as dtype objects). Names creates fields with labels.
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)
    return time, total, data
