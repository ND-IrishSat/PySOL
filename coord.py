''' coord.py

    Module for handling coordinate transformations

    Juwan Jeremy Jacobe
    University of Notre Dame
'''

import numpy as np

omega_Earth = np.array([0, 0, 7.3e-5]) # angular velocity of Earth / moving system with respect to ECI. Units of rad/s

# Implement ECI (or GCRF) --> ECEF (and or ITRF) for velocity and acceleration
def v_inert2v_rot(r_inert: np.ndarray, v_inert: np.ndarray, w_ri: np.ndarray = omega_Earth, v_org: np.ndarray = np.array([0, 0, 0])):
    ''' Transform v_rot to v_inert where v_rot is in a rotating frame 

        Args:
            r_inert: position in inertial frame (ECI)
            v_inert: velocity in inertial frame (ECI)
            w_ri: angular velocity of moving frame (w_Earth)
            v_org: linear velocity of origin of moving frame, default of 0
            
        Returns:
            r_ECEF
            v_ECEF
    '''

def v_rot2vinert()

