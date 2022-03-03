## WMM Model 2020
This folder contains a module implementing the World Magnetic Model 2020 in Python as well as its corresponding csv file containing the Gauss coefficietns needed 
for performing calculations. The primary usage for this module is to be calculate the B-Field due to the Earth at a certain coordinates (given in geodesic latitude,
longitude, and height from the ellipsoid) and at a certain time given in decimal years

This current implementation requires:
- numpy
- pyshtools
