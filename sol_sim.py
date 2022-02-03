"""
    sol_sim.py

    



"""

from time import strftime
import spacecraft as sp 
import orb_tools as ot 
import models

import datetime as dt
import numpy as np
import scipy.integrate as sci
import matplotlib.pyplot as plt
import constants

class Simulation():

    def __init__(self, model_name = 'Two-Body', TIME = None):
        """
            __init__ - initialization of PySol Simulation

                Opt. Args:
                    model_name (str) :
                        Dynamical model name 
                        Current dynamical models:
                            ['Two-Body']

                    TIME (dt.datetime object) : 
                        UTC time to initialize simulation
                        if None, sim time set to current UTC

        """

        # Define Constants used in simulation
        self.R_E = constants.R_E()

        # Generate dynamical model object
        self.model_nm = model_name
        model = models.Orbital_Models(model_name)

        # State integration function
        self.state_func = model.get_state_func()

        # Set simulation clock
        if TIME == None:
            self.t0  = dt.datetime.utcnow()
            self.time = self.t0 
        else:
            self.t0 = TIME 
            self.time = TIME

        # Intitialize spacecraft and sim time lists
        self.scs = []
        self.times = []

        pass

    def create_sc(self, OE_array, verbose = False, color = 'firebrick', name = None):
        """
            create_sc - function takes OE_array and generates sc object in memeory

                Args:
                    OE_array (OE_array obj): OE_array object to initialize sc
        
                Returns:
                    sc (sc obj): spacecraft object
        """

        # Time of init is taken from current sim time
        TIME = self.time

        sc = sp.Spacecraft(OE_array, TIME, verbose = verbose, color = color, name = name)
        self.scs.append(sc)

        return sc

    def propogate(self, DT, resolution = 10, tol = [1e-7, 1e-4], integrator = 'RK45', 
        event_func = None):
        """
            propogate - propogates sc and simulation forward in time using dynamical model

            Args:
                DT (datetime.timedelta) : time of integration
                resolution (float) : output time spacing in seconds
                tol (1x2 list) : tolerance of RK integrator
                integrator (str) : num. integration method to use 
                    ['RK23', 'RK45', 'DOP853', 'Radau']
                event_func (function) : function to record events

            Return:
                None
        """

        dt_seconds = DT.seconds
        n_outputs = int(dt_seconds/resolution)

        # Propogate each sc and record output to sc object
        for i, sc in enumerate(self.scs):
            sc_props = self.propogate_func(sc, dt_seconds, n_outputs, tol, integrator, event_func)

            # Add new state vectors to sc time list
            states = sc_props[0]
            self.scs[i].set_states(states)
            
            # Add new datetime objects to sc time list
            t_s = sc_props[1]
            times = []
            for s in t_s:
                delta_t = dt.timedelta(seconds  = s)
                times.append(self.time + delta_t)
            self.scs[i].set_times(times)
        
        # Set the sim time to the last sc recorded time
        self.time = times[-1]
        self.times = times

    def propogate_func(self, sc, dt_sec, n_outputs, tol, integrator, event_func):
        """
            propogate_func:
                function for use in the simulation.propogate_parallel, 
                    simulation.propogate methods

                    Args:
                        sc (sc object): sc object for propogation
                        tau_f (float): final non-dim time
                        n_outputs (int): number of integration outputs
                                        (NOT # of integration steps!) 
                        tol (tuple): absolute and relative intergation tolerance, set to
                            [1e-12, 1e-10] for high precision
                        event_func: func to record events

                    Returns:
                        s_new (6xN): sc states at each t_eval
                        t_new (N): evaluated times

                    Opt. Returns:
                        y_hits: states at event triggers
                        t_hits: times at event triggers 
                         
        """

        # Set up integration
        state = sc.state_mat.S_[-1]
        t_eval = np.linspace(0, dt_sec, n_outputs)

        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
        # Numerical integration 
        sol = sci.solve_ivp(
            fun = self.state_func,  # Integrtion func
            y0 = state,             # Initial state
            t_span = [0, dt_sec],   # Init/Final int time
            method = integrator,    # Int Algorithm
            t_eval = t_eval,        # Times to output
            max_step = 10,          # Max time step [s] in int
            atol = tol[0],
            rtol = tol[1],
            #args = [eval_STM, p_event, epoch],
            events = event_func
        )

        # propogated states and time arrays
        s_new = sol.y.T
        t_new = sol.t.T

        return s_new, t_new


    def plot_orbit(self, D2 = False, tau_f = None, Earth = True, 
        lims = [8000, 8000, 8000], IJK = True):

        plt.rcParams.update({'font.sans-serif': 'Helvetica'})

        xlim, ylim, zlim = lims

        if D2:
            fig = plt.figure(figsize = [10, 10])
            ax = fig.add_subplot()
            ax.set_aspect('equal')

            if Earth:
                earth = self.__earth_2d()
                ax.add_artist(earth)
        else:
            
            fig = plt.figure(figsize = [8, 8])
            ax = fig.add_subplot(projection = '3d')
            ax.set_xlim(-xlim, xlim)
            ax.set_ylim(-ylim, ylim)
            ax.set_zlim(-zlim, zlim)
            ax.set_box_aspect([1, ylim/xlim, zlim/xlim])

            ax.set_xlabel('X [km]')
            ax.set_ylabel('Y [km]')
            ax.set_zlabel('Z [km]')

            if Earth:
                earth = self.__earth_3d()
                ax.plot_wireframe(earth[0], earth[1], earth[2], 
                    color = 'mediumblue', alpha = 0.1, zorder = 0)

            if IJK:
                I = xlim*np.array([[0, 0, 0], [1, 0, 0]])
                J = ylim*np.array([[0, 0, 0], [0, 1, 0]])
                K = zlim*np.array([[0, 0, 0], [0, 0, 1]])

                plt.plot(I[:, 0], I[:, 1], I[:, 2], color = 'black')
                plt.plot(J[:, 0], J[:, 1], J[:, 2], color = 'black')
                plt.plot(K[:, 0], K[:, 1], K[:, 2], color = 'black')

        tf = self.scs[0].state_mat.times[-1]
        t0 = self.scs[0].state_mat.times[0]
        
        ax.set_title('PySOL | Dynamics: {} | '.format(self.model_nm) + t0.strftime('%Y/%m/%d') +
            '\n' + t0.strftime('%H:%M:%S – ') + tf.strftime('%H:%M:%S UTC'))
        for sc in self.scs:
            X = sc.state_mat.X
            Y = sc.state_mat.Y
            Z = sc.state_mat.Z
            ax.plot(X, Y, Z, color = sc.color, zorder = 2)
            ax.scatter(X[-1], Y[-1], Z[-1], s = 100, 
                fc = 'black', ec = sc.color, marker = '.', zorder = 3, label = sc.name)


        ax.legend()

    def plot_RADEC(self,):

        fig = plt.figure()
        ax = fig.add_subplot()

        for sc in self.scs:
            RADEC = sc.state_mat.to_RADEC()
            ax.plot(RADEC[:, 0], RADEC[:, 1])

        ax.xlabel('Right Ascension [deg]')
        ax.ylabel('Declination [deg]')
        
    def plot_XYZ(self):

        fig, axs = plt.subplots(nrows = 3, ncols = 1, figsize = [12, 6], sharex = True)

        plt.rcParams.update({'font.sans-serif': 'Helvetica'})

        for sc in self.scs:
            labels = ['X [km]', 'Y [km]', 'Z [km]']
            states = sc.state_mat.S_
            times = sc.state_mat.times
            for i, ax in enumerate(axs):
                if i < 1:
                    axs[i].plot(times, states[:, i], color = sc.color, label = sc.name)
                else:
                    axs[i].plot(times, states[:, i], color = sc.color)

                if i > 1:
                    axs[i].set_xlabel('Time [UTC]')

                ax.minorticks_on()
                ax.tick_params('both', which = 'major', length = 10, direction = 'in')
                ax.tick_params('both', which = 'minor', length = 5, direction = 'in')
                ax.set_ylabel(labels[i])

        axs[0].legend(loc = 0)

        
    def __earth_2d(self):

        R_e = 6371 #km
        earth = plt.Circle((0, 0), R_e, color = 'mediumblue')

        return earth

    def __earth_3d(self):
        """ 
            __earth_3d:
            produces 3D wireplot state vector of Earth in the synodic frame

                Args:
                    r (1x6 array): position of Earth in ndim synodic frame

                Returns 
                    earth (1x6 array): Earth wireframe state vector
        """

        R_e = 6371 #km

        u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]
        x = (np.cos(u)*np.sin(v))*R_e
        y = (np.sin(u)*np.sin(v))*R_e
        z = (np.cos(v))*R_e

        earth = np.array([x, y, z])

        return earth 


if __name__ == '__main__':

    sim = Simulation()


    OE1 = ot.OE_array(f = 7.5, a = 8000, e = 0.00068, i = 51.64, Om = 300, w = 74)
    sim.create_sc(OE_array= OE1, verbose = True, color = 'green', name = 'IRISAT')

    # OE2 = ot.OE_array(f = 0, a = 20_000, e = 0.36, i = 46, Om = 4, w = 8)
    # sim.create_sc(OE_array= OE2, verbose = True, name = 'Bridget\'s Orbit', color = 'Pink')

    # OE3 = ot.OE_array(f = 0, a = 5000, e = 0.1, i = 69, Om = 42, w = 63)
    # sim.create_sc(OE_array= OE3, verbose = True, color = 'darkorange', name = 'Owen')

    DT = dt.timedelta(hours = 1)
    sim.propogate(DT, resolution = 1)

    sim.plot_orbit(lims = [6000, 6000, 6000])

    sim.plot_RADEC()

    #sim.plot_XYZ()

    plt.show()