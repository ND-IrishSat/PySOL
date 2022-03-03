''' Your one and only WMM model. Based on the report attached below. This program needs the
    file WMMcoef.csv, which contains the values for the Gauss coefficients used to compute the
    magnetic field components.

    This current version relies on the pyshtools library to calculate the needed Schmidt semi normalized associated
    Legendre functions, in which the library actually loads from a .pyd file so it's hard to actually control how this is done >:(
    Future implementations of this model will include a built-in script for directly calculating these legendre functions for more
    efficiency rather than loading this pretty wack library
    
    The runtime of this implementation of the WMM model is:
        ~0.03-0.04 s to define the model
        ~0.02 to make a magnetic field calculation
    
    by Juwan Jeremy Jacobe
    
    References:
    WMM 2020 Report: https://www.ngdc.noaa.gov/geomag/WMM/soft.shtml
    
    Needed libraries:
        numpy
'''

import numpy as np
import time
import pyshtools.legendre as legendre

class WMMCoefficientLoader():
    ''' Class to load and hold model coefficients from reference epoch the WMM geomagnetic model for nth degree  
    model
    
    Arguments for initializing:
        degree (int): degree of spherical harmonic coefficients to use. Ranges of [1, 12]
        
    '''
    
    def __init__(self, degree):  

        # Determine length of array of coefficients based on degree, where for given nth degree
        # there are m coefficients from m = 0 to m = n
        self.degree = degree
        self.length = 0
        
        for i in range(1, degree+1):
            self.length += i + 1
        
        # Initialize gauss coefficients as empty arrays
        self.g_ = np.zeros(self.length)
        self.h_ = np.zeros(self.length)
        self.gdot_ = np.zeros(self.length)
        self.hdot_ = np.zeros(self.length)
        
    def read_coefficients(self, file_name):
        ''' This reads the Gauss coefficients from the .csv file within this directory to the class by using
        numpy's loadtxt function
        
        Args:
            file_name (str): name of csv file to read coefficients from
        '''
        
        data = np.loadtxt('WMMcoef.csv')
        
        self.g_ = data[:, 0]
        self.h_ = data[:, 1]
        self.gdot_ = data[:, 2]
        self.hdot_ = data[:, 3]      
    
    def shape_coefficients(self):
        ''' Shapes the n degree coefficients with m ({0, m]) terms into an n dimensional tuple, where
        each element is m=[0,n] number of coefficients stored inside a numpy array
        
        '''
        
        # Shape each of the coefficients
        self.g_ = self.shaper(self.g_)
        self.gdot_ = self.shaper(self.gdot_)
        self.h_ = self.shaper(self.h_)
        self.hdot_ = self.shaper(self.hdot_)
        
    def shaper(self, coefficient):
        ''' Shapes a single coefficient into an n dimensional tuple, where each element is m = [0, n] number of coefficients
        stored inside a numpy array
        
        Arguments:
            coefficient (1d numpy array): array to be shaped into n dimensional tuple
            
        Out:
            shaped_coefficient (tuple): a tuple where each element is an np-array of the m coefficients for a given nth degree
        '''
        
        # Initial parameters
        shaped_coefficient = []
        n = 0
        index = 0
        
        # Keep shaping until all n degree coefficients are shaped
        while (n < self.degree):
            row = np.array([])
            
            # Take the m -> [0, n] coefficients as one np array
            for m in range(n+2):
                row = np.append(row, coefficient[index+m])
                
            # Don't append if list is empty
            if n == 0:
                shaped_coefficient = [row]
            # Else, append
            else:
                shaped_coefficient.append(row)
            
            # Incrementation
            n += 1
            index += n + 1
        
        return shaped_coefficient
                
            
class WMM():
    ''' Class to hold the WMM model and to calculate geomagnetic field components using
    the model equations
    
    The file_name you will be loading is 'WMMcoef.csv.' The method that you would use is calc_gcc_components(), which uses the other
    methods defined here. Treating the whole class a blackbox, the input is LLA coordinates, an array where the first element is the geodetic latitude,
    the second element longtitude, and the third element is height from ellipsoid (a form of approximating altitude?). The output of this model will be
    the X, Y, and Z magnetic field components in the ellipsoidal reference frame at that point 
    '''
    
    def __init__(self, degree, file_name):
        ''' Model initialization
        
        Args:
            degree (int): the degree to which you want to calculate the model for
            file_name (str): the path + file name you are reading Gauss coefficient values from
        '''
        
        # Load Gauss coefficients from file
        self.coef = WMMCoefficientLoader(degree)
        self.coef.read_coefficients(file_name)
        self.coef.shape_coefficients()
        self.degree = degree
        self.degree_max = 12
        
        ''' NOTE:
        As self.coef is a WMMCoefficientLoader class, the respective Gauss coefficients at the epoch time
        can be accessed via:
            self.coef.g_
            self.coef.gdot_
            self.coef.h_
            self.coef.hdot_
        '''
        
        # Other coefficients, you may want to use your custom coefficient class or whatevs
        self.a = 6371200 # geomagnetic reference radius, m
        self.t_0 = 2020.0 # epoch reference time, in decimal years
        self.A = 6378137 # semi-major axis, m
        self.reciprocal_f = 298.257223563 # 1/f, reciprocal of f
        self.f = 1.0/self.reciprocal_f
        self.e_earth = np.sqrt(self.f*(2-self.f))
            
    def calc_LLA2GCC(self, LLA, t, degrees=False):
        ''' Method to convert inputted geodetic latitude, longtitude, altitude into geocentric spherical coordinates,
        which will be stored inside the self of the model
        
        Args:
            LLA (numpy.array): A numpy array where the first element is the geodetic latitude, the second
            element is the longitude, and the third is the height above the ellipsoid. Height/altitude should be
            in meters
            
            t (number): time, given in decimal years
            
        The geocentric spherical coordinates is saved as an np array, where the first element is the goecentric latitude,
        the second is longitude, and the third is geocentric radius. The geodetic coorindates are also saved into
        the model for calculation purposes. The latitude and longtitude are stored 
        
        NOTE: DO NOT USE THIS YOURSELF
        '''
        
        # Store LLA vector to temporary variables
        lat_gd = LLA[0] 
        lon = LLA[1] # Not converted to radians as not used in calculations
        h_ellp = LLA[2]
        
        if degrees:
            lat_gd *= np.pi/180.0
            lon *= np.pi/180.0
            
        # Calculate radius of curvature of the prime vertical for given coordinates
        self.R_c = self.A / (np.sqrt(1-(self.f*(2-self.f))*(np.sin(lat_gd)**2)))
        
        # Radius in x and y cartesian geocentric (ECEF) plane
        p = (self.R_c + h_ellp)*np.cos(lat_gd)
        
        # Calculate z cartesian geocentric (ECEF) coordinate
        z = (self.R_c*(1-self.e_earth**2) + h_ellp)*np.sin(lat_gd)
        
        # Get geocentric radius from p and z
        r = np.sqrt(p**2 + z**2)
        
        # Get geocentric latitude
        lat_gc = np.arcsin(z/r)
        
        
        # Store geodesic vector, geodetic vector, as well as corresponding time stamp
        self.GDC = np.array([lat_gd, lon, h_ellp])
        self.GCC = np.array([lat_gc, lon, r])
        self.t = t
        
    def determine_coefficients(self):
        '''Determine gauss coefficients for desired time using equations:
        
        g^m_n(t) = g^m_n(t_0) + (t-t_0) g_dot^m_n(t_0)
        h^m_n(t) = h^m_n(t_0) + (t-t_0) h_dot^m_n(t_0)
        
        where g^m_n(t_0) and h^m_n(t_0) are the coeffients for the reference epoch, t_0 = 2020.0, of which
        the values are stored in the file and in the WMMCoefficientLoader class. The quantities g_dot^m_n(t_0) and
        h_dot^m_n(t_0) are the secular variation coefficients.
        
        This function also calculates the array of Schmidt semi-normalized associated Legrende functions of nth
        order and mth order, where P_m,n is stored in the ( n*(n+1)/2 + m )th index of the matrix 
        
        --------------------------------------------------------------------------------------------------
        IMPORTANT TO NOTE:
        The way this is being done right now is quite honestly inefficient, this will be improved upon by flattening it to
        a 1D numpy array and having a smart iterator
        
        '''
        
        # Empty lists
        self.g_t = []
        self.h_t = []
        self.g_dot_t = []
        self.h_dot_t = []
        
        # Calculate Gauss coefficients for desired time using model coefficients set from epoch 
        for i in range(self.degree):
            
            # Set the m number of coefficients for a given order n
            g_coefs = np.array(self.coef.g_[i])
            h_coefs = np.array(self.coef.h_[i])
            
            # Pull the rate of change for the m coefficients for a given order n
            g_coefs_dot = np.array(self.coef.gdot_[i])
            h_coefs_dot = np.array(self.coef.hdot_[i])
            
            # Take into account time dependency based on difference between given time and epoch time of the model, which is 2020.0 in decimal years
            g_coefs += (self.t - self.t_0)*g_coefs_dot
            h_coefs += (self.t - self.t_0)*h_coefs_dot
            
            # Append to the list 
            self.g_t.append(g_coefs)
            self.h_t.append(h_coefs)
            self.g_dot_t.append(g_coefs_dot)
            self.h_dot_t.append(h_coefs_dot)
            
        # Calculate the array of Schmidt semi normalized associated legendre functions for given latitude and given order    
        self.legendre = legendre.PlmSchmidt(self.degree+1, np.sin(self.GCC[0]))
        
                  
    def calc_gcc_components(self, LLA, t, degrees=False):
        ''' Determine the field vector components X', Y', and Z' in geocentric coordinates, save to model
        and also output:
        
        Args:
            LLA (np.array): geodetic coordinates, first element is geodetic latitude, second element is longitude,
            and third element is height above ellipsoid
            t (number): Associated time of coordinates in decimal years, where t = 2020.0 refers to the year 2020
            degrees (bool): Whether or not inputted LLA coordinates are in degrees or radians. True for degrees
            
        Out:
            B_gcc (np.array): geocentric magnetic field vector [X', Y', Z']
        '''
        
        # Convert inputted geodetic latitude, longtitude, altitude coordinates into spherical geocentric coordinates
        # of geocentric latitude, longtitude, and radius from center
        self.calc_LLA2GCC(LLA, t, degrees)
        self.determine_coefficients()
        
        # Intermediate values, pulling lon, lat, r, and other values from self for better readability
        lat = self.GCC[0]     # radians
        lon = self.GCC[1]     # radians
        r = self.GCC[2]       # m
        a = self.a            # m
        lat_gd = self.GDC[0]  # radians
        h = self.GDC[2]       # m
        
        # Field vector components
        B_x = 0               # nT
        B_y = 0               # nT
        B_z = 0               # nT
        
        # Equations 10-12 of the WMM 2020 Report, the sum over n = 1 to n = 12
        for n in range(1, self.degree+1):
        
            # Initialize sum term over n as 0
            B_x_n = 0
            B_y_n = 0
            B_z_n = 0
            
            # The sum over m from m = 0 to m = n
            for m in range(n+1):
            
                # Hold Gaussian coefficients in intermediate values for convenience and to not access arrays
                # like 20 times :)
                g_t = self.g_t[n-1][m] # We do n-1, but m because n is counted from 1 but m is counted from 0 for math purposes
                h_t = self.h_t[n-1][m]
                g_t_dot = self.g_dot_t[n-1][m]
                h_t_dot = self.h_dot_t[n-1][m]
               
                # SSNA Legrende function of sin(latitude) with degree being n and order being m
                Lp_mn = self.legendre[int(n*(n+1)/2 + m)]                

                # SSNA Legrende function of sin(latitude) with degree being n+1 and order being m
                Lp_mn1 = self.legendre[int((n+1)*(n+2)/2 + m)]
                
                # Derivative of SSNA Legrende function with respect to latitude
                Lp_derivative =  ((n+1)*np.tan(lat) * Lp_mn) - np.sqrt((n+1)**2 - m**2) * (1/np.cos(lat)) * Lp_mn1

                # Sump up over m
                B_x_n += (g_t * np.cos(m*lon) + h_t*np.sin(m*lon)) * Lp_derivative
                B_y_n += m*(g_t * np.sin(m*lon) - h_t*np.cos(m*lon)) * Lp_mn
                B_z_n += (g_t * np.cos(m*lon) + h_t*np.sin(m*lon)) * Lp_mn
                
                
            # Multiply with factor outside sum over m
            B_x_n *= (a/r)**(n+2)
            B_y_n *= (a/r)**(n+2)
            B_z_n *= (n+1)*(a/r)**(n+2)
            
            
            # Add the nth instance to the total component
            B_x += B_x_n
            B_y += B_y_n
            B_z += B_z_n
        
        # Final factors, outside of summation
        B_x *= -1
        B_y *= (1/np.cos(lat))
        B_z *= -1
        
        # Save geocentric magnetic field 
        self.GCBfield = np.array([B_x, B_y, B_z])
        
        # Rotate to ellipsoidal reference frame
        B_x_ellip = B_x * np.cos(lat - lat_gd) - B_z * np.sin(lat - lat_gd)
        B_y_ellip = B_y
        B_z_ellip = B_x * np.sin(lat - lat_gd) + B_z * np.cos(lat - lat_gd)
        
        # Save B field in elliptical reference frame
        self.Bfield_ellip = np.array([B_x_ellip, B_y_ellip, B_z_ellip])
        
    def get_Bfield(self):
        ''' Function to return the B field components in the ellipsoidal reference frame '''
        return self.Bfield_ellip
       
# Testing WMM model
if __name__ == "__main__":

    time1 = time.time()
    ## Test case of model
    n = 12
    WMM_model = WMM(n, 'WMMcoef1.csv')
    time2 = time.time()
    print(f'Time to create model: {time2-time1}\n')
    
    # Test example
    lat = -80.0 # latitude
    lon = 240.0 # longitude
    h = 100000.0 # altitude in m
    t = 2022.5 # year
    
    # Test calculation
    print(f'Computing for latitude = {lat}, longitude = {lon}, altitude = {h} m: ')
    time1 = time.time()
    WMM_model.calc_gcc_components(np.array([lat, lon, h]), t, degrees=True)
    time2 = time.time()
    print(f'Time to compute: {time2-time1}')
    print(WMM_model.get_Bfield())