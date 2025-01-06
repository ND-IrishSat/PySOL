''' 
    Module containing analytical models for the B field of the frames of the Helmholtz cage + Helmholtz cage.
    
    by Juwan Jeremy Jacobe
    Refactored by Andres Perez
    University of Notre Dame
    IrishSat OAT Lab -> IrishSat Goat Lab
'''

import numpy as np
import matplotlib.pyplot as plt

#####################################################################
#      Global Frame Classes                                         #
#####################################################################
# The original multi-frame classes have been refactored into a global frame class. The equations
# for the B-field are derived from the Biot-Savart Law.
#
# Currently, it is implemented as one major class. However, future iterations could explore the idea
# of having a Wire class to make up frames, as the equations are mostly the same!
#
# While this may not provide significant functional optimization, it offers a valuable learning opportunity
# to understand the Helmholtz cage, Helmholtz coils, and magnetic fields in general.
#
# NOTE: The B-field calculations from the frames have been independently confirmed to be correct for
# simple cases, such as the B-field at the center of the frame.

class GlobalFrame():
    def __init__(self, length, num, displacement, frame):
        
        # Checks if frame is initialized correctly
        if frame not in ["x", "y", "z"]:
            raise ValueError("frame must be 'x', 'y', or 'z'")

        # Constants:
        self.length = length
        self.N = num
        self.D = displacement 
        self.frame = frame
        
        # Constant in front of Biot-Savart Law: mu_0*N/(4*pi) and Permeability Constant
            # We're leaving out current (A), I, for modularity purposes
        self.mu_0 = 1.256637062e-6 
        self.const = self.mu_0*num/(4*np.pi)
    
        # Set displacements based on frame orientation
        if frame == "x":
            # Frame parallel to z-y plane
            # x displacement the frame from the origin
            self.primary_disp = displacement
            self.side1_disp = length / 2    # y displacement of wire 1
            self.side2_disp = length / 2    # z displacement of wire 2
            self.side3_disp = -length / 2   # y displacement of wire 3
            self.side4_disp = -length / 2   # z displacement of wire 4
        elif frame == "y":
            # Frame parallel to x-z plane
            self.primary_disp = displacement
            self.side1_disp = -length / 2   # x displacement of wire 1
            self.side2_disp = length / 2    # z displacement of wire 2
            self.side3_disp = length / 2    # x displacement of wire 3
            self.side4_disp = -length / 2   # z displacement of wire 4
        elif frame == "z":
            # Frame parallel to x-y plane
            self.primary_disp = displacement
            self.side1_disp = length / 2    # y displacement of wire 1
            self.side3_disp = -length / 2   # y displacement of wire 3
            self.side2_disp = -length / 2   # x displacement of wire 2
            self.side4_disp = length / 2    # x displacement of wire 4
    
    def get_B_factors(self,x,y,z):
        if self.frame == "x":
            # integrated over z'
            B1_integral_term1 = (self.length - 2 * z) / np.sqrt(4 * ((x - self.primary_disp) ** 2 + (y - self.side1_disp) ** 2) + (self.length - 2 * z) ** 2)
            B1_integral_term2 = (self.length + 2 * z) / np.sqrt(4 * ((x - self.primary_disp) ** 2 + (y - self.side1_disp) ** 2) + (self.length + 2 * z) ** 2)
            B1_integral = (B1_integral_term1 + B1_integral_term2) / ((x - self.primary_disp) ** 2 + (y - self.side1_disp) ** 2)
            
            # integrated over y'
            B2_integral_term1 = (self.length - 2 * y) / np.sqrt(4 * ((x - self.primary_disp) ** 2 + (z - self.side2_disp) ** 2) + (self.length - 2 * y) ** 2)
            B2_integral_term2 = (self.length + 2 * y) / np.sqrt(4 * ((x - self.primary_disp) ** 2 + (z - self.side2_disp) ** 2) + (self.length + 2 * y) ** 2)
            B2_integral = -(B2_integral_term1 + B2_integral_term2) / ((x - self.primary_disp) ** 2 + (z - self.side2_disp) ** 2)
            
            # integrated over z'
            B3_integral_term1 = (-self.length + 2 * z) / np.sqrt(4 * ((x - self.primary_disp) ** 2 + (y - self.side3_disp) ** 2) + (self.length - 2 * z) ** 2)
            B3_integral_term2 = (-self.length - 2 * z) / np.sqrt(4 * ((x - self.primary_disp) ** 2 + (y - self.side3_disp) ** 2) + (self.length + 2 * z) ** 2)
            B3_integral = (B3_integral_term1 + B3_integral_term2) / ((x - self.primary_disp) ** 2 + (y - self.side3_disp) ** 2)
            
            # integrated over y'
            B4_integral_term1 = (-self.length + 2 * y) / np.sqrt(4 * ((x - self.primary_disp) ** 2 + (z - self.side4_disp) ** 2) + (self.length - 2 * y) ** 2)
            B4_integral_term2 = (-self.length - 2 * y) / np.sqrt(4 * ((x - self.primary_disp) ** 2 + (z - self.side4_disp) ** 2) + (self.length + 2 * y) ** 2)
            B4_integral = -(B4_integral_term1 + B4_integral_term2) / ((x - self.primary_disp) ** 2 + (z - self.side4_disp) ** 2)

        elif self.frame == "y":
            # integrated over z'
            B1_integral_term1 = (self.length - 2 * z) / np.sqrt(4 * ((y - self.primary_disp) ** 2 + (x - self.side1_disp) ** 2) + (self.length - 2 * z) ** 2)
            B1_integral_term2 = (self.length + 2 * z) / np.sqrt(4 * ((y - self.primary_disp) ** 2 + (x - self.side1_disp) ** 2) + (self.length + 2 * z) ** 2)
            B1_integral = (B1_integral_term1 + B1_integral_term2) / ((y - self.primary_disp) ** 2 + (x - self.side1_disp) ** 2)
            
            # integrated over x'
            B2_integral_term1 = (self.length - 2 * x) / np.sqrt(4 * ((y - self.primary_disp) ** 2 + (z - self.side2_disp) ** 2) + (self.length - 2 * x) ** 2)
            B2_integral_term2 = (self.length + 2 * x) / np.sqrt(4 * ((y - self.primary_disp) ** 2 + (z - self.side2_disp) ** 2) + (self.length + 2 * x) ** 2)
            B2_integral = (B2_integral_term1 + B2_integral_term2) / ((y - self.primary_disp) ** 2 + (z - self.side2_disp) ** 2)
            
            # integrated over z'
            B3_integral_term1 = (-self.length + 2 * z) / np.sqrt(4 * ((y - self.primary_disp) ** 2 + (x - self.side3_disp) ** 2) + (self.length - 2 * z) ** 2)
            B3_integral_term2 = (-self.length - 2 * z) / np.sqrt(4 * ((y - self.primary_disp) ** 2 + (x - self.side3_disp) ** 2) + (self.length + 2 * z) ** 2)
            B3_integral = (B3_integral_term1 + B3_integral_term2) / ((y - self.primary_disp) ** 2 + (x - self.side3_disp) ** 2)

            # integrated over x'
            B4_integral_term1 = (-self.length + 2 * x) / np.sqrt(4 * ((y - self.primary_disp) ** 2 + (z - self.side4_disp) ** 2) + (self.length - 2 * x) ** 2)
            B4_integral_term2 = (-self.length - 2 * x) / np.sqrt(4 * ((y - self.primary_disp) ** 2 + (z - self.side4_disp) ** 2) + (self.length + 2 * x) ** 2)
            B4_integral = (B4_integral_term1 + B4_integral_term2) / ((y - self.primary_disp) ** 2 + (z - self.side4_disp) ** 2)

        elif self.frame == "z":
            # integrated over x'
            B1_integral_term1 = (-self.length + 2 * x) / np.sqrt(4 * ((y - self.side1_disp) ** 2 + (z - self.primary_disp) ** 2) + (self.length - 2 * x) ** 2)
            B1_integral_term2 = (-self.length - 2 * x) / np.sqrt(4 * ((y - self.side1_disp) ** 2 + (z - self.primary_disp) ** 2) + (self.length + 2 * x) ** 2)
            B1_integral = (B1_integral_term1 + B1_integral_term2) / ((y - self.side1_disp) ** 2 + (z - self.primary_disp) ** 2)
            
            # integrated over y'
            B2_integral_term1 = (self.length - 2 * y) / np.sqrt(4 * ((x - self.side2_disp) ** 2 + (z - self.primary_disp) ** 2) + (self.length - 2 * y) ** 2)
            B2_integral_term2 = (self.length + 2 * y) / np.sqrt(4 * ((x - self.side2_disp) ** 2 + (z - self.primary_disp) ** 2) + (self.length + 2 * y) ** 2)
            B2_integral = -(B2_integral_term1 + B2_integral_term2) / ((x - self.side2_disp) ** 2 + (z - self.primary_disp) ** 2)
            
            # integrated over x'
            B3_integral_term1 = (self.length - 2 * x) / np.sqrt(4 * ((y - self.side3_disp) ** 2 + (z - self.primary_disp) ** 2) + (self.length - 2 * x) ** 2)
            B3_integral_term2 = (self.length + 2 * x) / np.sqrt(4 * ((y - self.side3_disp) ** 2 + (z - self.primary_disp) ** 2) + (self.length + 2 * x) ** 2)
            B3_integral = (B3_integral_term1 + B3_integral_term2) / ((y - self.side3_disp) ** 2 + (z - self.primary_disp) ** 2)
            
            # integrated over y'
            B4_integral_term1 = (-self.length + 2 * y) / np.sqrt(4 * ((x - self.side4_disp) ** 2 + (z - self.primary_disp) ** 2) + (self.length - 2 * y) ** 2)
            B4_integral_term2 = (-self.length - 2 * y) / np.sqrt(4 * ((x - self.side4_disp) ** 2 + (z - self.primary_disp) ** 2) + (self.length + 2 * y) ** 2)
            B4_integral = -(B4_integral_term1 + B4_integral_term2) / ((x - self.side4_disp) ** 2 + (z - self.primary_disp) ** 2)
            
        self.B1_factor = self.const * B1_integral 
        self.B2_factor = self.const * B2_integral 
        self.B3_factor = self.const * B3_integral 
        self.B4_factor = self.const * B4_integral
    
    def get_Bfield(self, I, r):
        ''' 
            get_BfieldMM() - Takes the wire's magnetic field generated in 'get_B_factors()'
                            and multiplies it by the current. This calculates the coil's x, y, and z 
                            magnetic field components.
            
            Args:
                I (A)   : current
                r       : a position vector represented by a tuple/vector/array holding the values of x, y, and z,
        '''
        # r = [x, y, z]
        x, y, z = r
        self.get_B_factors(x, y, z)

        if self.frame == "x":
            # Frame parallel to the y-z plane
            B_x = I * (-self.B1_factor * (y - self.side1_disp) + self.B2_factor * (z - self.side2_disp) +
                    -self.B3_factor * (y - self.side3_disp) + self.B4_factor * (z - self.side4_disp))
            B_y = I * (self.B1_factor * (x - self.primary_disp) + self.B3_factor * (x - self.primary_disp))
            B_z = I * (-self.B2_factor * (x - self.primary_disp) + -self.B4_factor * (x - self.primary_disp))

        elif self.frame == "y":
            # Frame parallel to the x-z plane
            B_x = I * (-self.B1_factor * (z - self.side2_disp) + -self.B3_factor * (z - self.side4_disp))
            B_y = I * (self.B1_factor * (x - self.side1_disp) + -self.B2_factor * (z - self.side2_disp) +
                    self.B3_factor * (x - self.side3_disp) + -self.B4_factor * (z - self.side4_disp))
            B_z = I * (self.B2_factor * (y - self.primary_disp) + self.B4_factor * (y - self.primary_disp))

        elif self.frame == "z":
            # Frame parallel to the x-y plane
            B_x = I * (-self.B2_factor * (z - self.primary_disp) + -self.B4_factor * (z - self.primary_disp))
            B_y = I * (self.B1_factor * (z - self.primary_disp) + self.B3_factor * (z - self.primary_disp))
            B_z = I * (self.B1_factor * (y - self.side1_disp) + -self.B2_factor * (x - self.side2_disp) +
                    self.B3_factor * (y - self.side3_disp) + -self.B4_factor * (x - self.side4_disp))

        return np.array([B_x, B_y, B_z])

###############################################################
#     Helmholtz Cage Analytical Model                         #
###############################################################
#  This class allows one to solve for the B-field due to a Helmholtz cage at position r
#  where the origin lies in the very center of this Helmholtz cage :>

class HelmholtzCage():
    def __init__(self, length, num, x_disp, y_disp, z_disp):
        ''' Initialize the cage object. Assumes that each frame is of equal side length + equal num of coils wrapped
        around the frame + the x, y, z frames are displaced equally from the origin.
            
            Args:
                length: the length of the frames making up the cage
                num: the number of coils wrapped around the frames
                x_disp: displacements of frames generating B_x from origin, where the origin is the center of the cage 
                y_disp: displacements of frames generating B_y from origin, where the origin is the center of the cage
                z_disp: displacements of frames generating B_z from origin, where the origin is the center of the cage
        '''
        
        # Frame generating +B_x, displaced -x_disp from origin)
        self.frame_xminus = GlobalFrame(length=length, num=num, displacement=-x_disp, frame='x')

        # Frame generating -B_x, displace +x_disp from origin
        self.frame_xplus = GlobalFrame(length=length, num=num, displacement=x_disp, frame='x')

        # Frame generating +B_y, displaced -y_disp from origin
        self.frame_yminus = GlobalFrame(length=length, num=num, displacement=-y_disp, frame='y')
        
        # Frame generating -B_y, displaced +y_disp from origin
        self.frame_yplus = GlobalFrame(length=length, num=num, displacement=y_disp, frame='y')
        
        # Frame generating +B_z, displaced -z_disp from origin
        self.frame_zminus = GlobalFrame(length=length, num=num, displacement=-z_disp, frame='z')
        
        # Frame generating -B_z, displaced +z_disp from origin
        self.frame_zplus = GlobalFrame(length=length, num=num, displacement=z_disp, frame='z')
        
        # Create tuple of all frames together, for coding elegance purposes :>
        self.cage = (self.frame_xminus, self.frame_xplus, self.frame_yminus, self.frame_yplus, self.frame_zminus, self.frame_zplus)
        
    def calc_Bfield(self, r, I_xminus=0.1, I_xplus= 0.1, I_yminus= 0.1, I_yplus= 0.1, I_zminus=0.1, I_zplus= 0.1):
        ''' Calculate B field due to Helmholtz cage at position r
        
            Args:
                r (np.array): array of x, y, z positions at which to calculate B field for
                I_xminus (number): current running through Frame -X   (A)
                I_xplus (number): current running through Frame +X    (A) 
                I_yminus (number): current running through Frame -Y   (A)
                I_yplus (number): current running through Frame +Y    (A)
                I_zminus (number): current running through Frame -Z   (A)
                I_zplus (number): current running through Frame +Z    (A)
            
            Out:
                B_field (np.array): array of magnetic field components, Bx, By, Bz
        '''
        
        B_field = np.zeros(3)
        
        currents = (I_xminus, I_xplus, I_yminus, I_yplus, I_zminus, I_zplus)
        
        for i, frame in enumerate(self.cage):
            B_field = B_field + frame.get_Bfield(currents[i], r)
            
        return B_field
        
    def calc_frame_vertices(self, frame):
        ''' Find the four vertices of the frame:
        
            Args:
                frame: a FrameX, FrameY, or FrameZ object
        '''
        if (frame==self.frame_xminus or frame==self.frame_xplus):
            vertix1 = [frame.primary_disp, frame.side1_disp, frame.side2_disp]
            vertix2 = [frame.primary_disp, frame.side3_disp, frame.side2_disp]
            vertix3 = [frame.primary_disp, frame.side3_disp, frame.side4_disp]
            vertix4 = [frame.primary_disp, frame.side1_disp, frame.side4_disp]
        if (frame==self.frame_yminus or frame==self.frame_yplus):
            vertix1 = [frame.side1_disp, frame.primary_disp, frame.side2_disp]
            vertix2 = [frame.side3_disp, frame.primary_disp, frame.side2_disp]
            vertix3 = [frame.side3_disp, frame.primary_disp, frame.side4_disp]
            vertix4 = [frame.side1_disp, frame.primary_disp, frame.side4_disp]
        if (frame==self.frame_zminus or frame==self.frame_zplus):
            vertix1 = [frame.side2_disp, frame.side1_disp, frame.primary_disp]
            vertix2 = [frame.side2_disp, frame.side3_disp, frame.primary_disp]
            vertix3 = [frame.side4_disp, frame.side3_disp, frame.primary_disp]
            vertix4 = [frame.side4_disp, frame.side1_disp, frame.primary_disp]
                
        return np.array([vertix1, vertix2, vertix3, vertix4, vertix1])
        
    def draw_cage(self, ax):
        ''' Draws cage on a figure
            Arg: 
                ax (Figure): figure object to plot cage onto
        '''
        for frame in self.cage:
            verts = self.calc_frame_vertices(frame)
            print(verts)
            ax.plot(verts[:, 0], verts[:, 1], verts[:, 2], marker = 'o', color = 'blue')
            
    def plot_Bfield(self, xlim, ylim, zlim, grid_num = 4, I_xminus=0.1, I_xplus= 0.1, I_yminus=0, I_yplus= 0, I_zminus= 0, I_zplus= 0):
        ''' Plot a 3D vector field of the B field due to the Helmholtz cage
        
            Args:
                xlim (number): + and - x limits to make grid for  (m)
                ylim (number): + and - y limits to make grid for  (m)    
                zlim (number): + and - z limits to make grid for  (m)
                grid_num (int): number of points per axis in mesh grid
                I_xminus (number): current running through Frame -X   (A)
                I_xplus (number): current running through Frame +X    (A) 
                I_yminus (number): current running through Frame -Y   (A)
                I_yplus (number): current running through Frame +Y    (A)
                I_zminus (number): current running through Frame -Z   (A)
                I_zplus (number): current running through Frame +Z    (A)
        '''
        
        # Create 3D mesh grid for x, y, z over limits
        x = np.linspace(-xlim, xlim, grid_num)
        y = np.linspace(-ylim, ylim, grid_num)
        z = np.linspace(-zlim, zlim, grid_num)
        x_mesh, y_mesh, z_mesh = np.meshgrid(x, y, z)
        
        # Initialize B field components as a 3D mesh
        B_cage_x = np.zeros((grid_num, grid_num, grid_num))
        B_cage_y = np.zeros((grid_num, grid_num, grid_num))
        B_cage_z = np.zeros((grid_num, grid_num, grid_num))
        
        # Find Bx, By, Bz for every x, y, z combination
        # Is what I'm doing naughty here? Yes. Is there another way to do it? Probably. Can it be easily implexmented? Who knows?
        for i in range(grid_num):
            for j in range(grid_num):
                for k in range(grid_num):
                    B_field = cage.calc_Bfield(np.array([x[i], y[j], z[k]]), I_xminus=I_xminus, I_xplus= I_xplus, I_yminus= I_yminus, I_yplus= I_yplus, I_zminus= I_zminus, I_zplus= I_zplus)
                    B_cage_x[i, j, k] = B_field[0]
                    B_cage_y[i, j, k] = B_field[1]
                    B_cage_z[i, j, k] = B_field[2]
        
        # Draw vector field
        ax = plt.figure().add_subplot(projection = '3d')
        ax.quiver(x_mesh, y_mesh, z_mesh, B_cage_x, B_cage_y, B_cage_z, length = 0.1, normalize = True, color = 'black')
        
        # Draw cage
        self.draw_cage(ax)
        
        # Labels
        ax.set_title('B-Field of the Helmholtz Cage')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_zlabel('z (m)')
        plt.show()
        
if __name__ == "__main__":
    # Test case
    L = 1
    num = 80

    # Set displacements of frame from origin (center of cage)
    x_disp = 0.25
    y_disp = 0.25
    z_disp = 0.25
    I = 1
    r = np.array([0.0, 0.0, 0.0])
    lim = 0.4
    
    cage = HelmholtzCage(L, num, x_disp, y_disp, z_disp)
    cage.plot_Bfield(lim, lim, lim)
    
