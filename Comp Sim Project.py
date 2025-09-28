import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import json
plt.rcParams.update({'font.size': 26})

class Body:
    
    def __init__(self, name ,mass, rx, vy, colour, display_radius):
        """
        Initialise a celestial body with given properties.

        Parameters
        ----------
        name : str
            Name of the celestial body.
        mass : float
            Mass of the body in Earth masses.
        rx : float
            Initial x-coordinate position.
        vy : float
            Initial velocity in the y-direction.
        colour : str
            Display color of the body.
        display_radius : float
            Radius used for visual representation (not physical size).
        neg_quad : boolian 
            States whether planet is in the bottom right quadrant of the coordinate system relative to the Sun

        """
        self.name = name
        self.mass = mass
        self.r = np.array([rx, 0.0])
        self.vy = vy
        self.v = np.array([0.0, self.vy])
        self.a = np.array([0.0, 0.0])
        self.a_prev = None
        self.colour = colour
        self.display_radius = display_radius
        self.neg_quad = False
        self.period_val = None
    
    def accelerate(self, other_bodies, G):
        """
        Calculate the acceleration of the body due to gravitational attraction 
        from other celestial bodies.

        """
        self.a = np.array([0.0, 0.0])
        for body in other_bodies:
            if self.name != body.name: #avoids the body being accelerated due to itself
                r_ba = self.r - body.r #relative position vecotr between two bodies
                a = -G * body.mass * r_ba / norm(r_ba)**3
                self.a += a
        return self.a

    def GPE(self, other_bodies, G):
        """
        Compute the gravitational potential energy of the body due to 
        interactions with other celestial bodies.

        """
        self.grav_energy = 0
        for body in other_bodies:
            if self.name != body.name: #avoids calculating the GPE of a body due to itself
                r_ba = self.r - body.r #relative position vector between two bodies 
                self.grav_energy += G * self.mass * body.mass / norm(r_ba)
        return self.grav_energy
                
    def KineticEnergy(self):
        """
        Compute the kinetic energy of the body.

        """
        self.KE = 1/2 * self.mass * (norm(self.v))**2
        return self.KE
    
    def period(self, sun, time):
        """
        Estimate the orbital period of the body based on its position 
        relative to the Sun.

        Notes
        -----
        This function detects when the body is in the lower right quadrant of a coordinate system with the 
        Sun at the origin before calculating the difference in the time between positive x-axis crossings 
        and therefore approximating its orbital period. It prints the calculated period.

        """
        r_sun = sun.r
        r_planet = self.r
        r_rel = r_planet - r_sun #calculate postiion of planet relative to the sun becuase sun drifts
        if self.name != sun.name: #this is to ensure the period of the sun is not calculated
            if r_rel[0] > 0 and r_rel[1] < 0: #detects if x coordinate is positive and y coordinate is negative
                self.neg_quad = True
            if self.neg_quad == True and r_rel[1] > 0 and not self.period_val:
                #calculates period if the body's y coordinate changes from negtive to positive after being in the lower right quadrant and if there is not already a value for the period. 
                self.period_val = time
                print(f"The orbit of {self.name} is: {time}")
                

class Satellite:
    def __init__(self, position, vy):
        """
        Initialise a satellite with an initial position and velocity.

        """
        self.r = position
        self.vy = vy
        self.v = np.array([0.0, self.vy])
        self.a = np.array([0.0, 0.0])
        self.a_prev = None #set to None so that Beeman integration method can't use its value before it has one
        self.closest_approaches = [np.inf, np.inf]
        
    def accelerate(self, other_bodies, G):
        """
        Calculate the gravitational acceleration of the satellite due to other 
        celestial bodies.

        """
        self.a = np.array([0.0, 0.0])
        for body in other_bodies:
            r_ba = self.r - body.r #relative position vector between two bodies
            a = -G * body.mass * r_ba / norm(r_ba)**3 
            self.a += a
        return self.a
    
    def closestDistance(self, bodies, time):
        """
        Track the closest approach of the satellite to target celestial bodies.

        Notes
        -----
        If the satellite gets within 0.014 AU of a body after 0.5 years, it prints the 
        time taken and the closest distance reached.

        """
        for i in range(len(self.closest_approaches)): #loops through list of two so calculation happens for both Mars and Earth
            distance = bodies[i].r - self.r
            
            #if the distance is less than the stored value of closest approach and more than 0.5 years has passed, change value of closest approach to distnace calculated
            if norm(distance) < norm(self.closest_approaches[i]) and time > 0.5: 
                self.closest_approaches[i] = distance
                if norm(self.closest_approaches[i]) < 0.014:
                    print(f"Time taken from launch to reach {bodies[i].name} is: {time}years")
                    print(f"The closest approach to {bodies[i].name} is: {norm(self.closest_approaches[i])}au ")
        
class Simulation:
    """
    A class to simulate the motion of celestial bodies and satellites 
    under the influence of gravity using numerical integration.
    
    """
    
    def __init__(self):
        """
        Initialises the simulation by loading parameters from a JSON file, 
        setting up simulation parameters, and creating a list of celestial bodies.

        """
        with open('parameters_solar.json') as f:
            parameters_solar = json.load(f)
        self.no_steps = parameters_solar['num_iterations']
        self.dt = parameters_solar['timestep']
        self.G = parameters_solar['grav_const']
        self.bodies = []
        self.total_energy = []
        self.satellite_xpositions = []
        self.satellite_ypositions = []
        self.satellites = None #set to none initially so it is easy to ignore when not being used
        self.frame_number = 0
        for body in parameters_solar['bodies']:
            if body['name'] == "Sun":
                self.bodies.append(Body(
                    body['name'], 
                    body['mass'],
                    body['orbital_radius'],
                    0,
                    body['colour'],
                    body['display_radius']))
                m_1 = self.bodies[0].mass
                r_1 = self.bodies[0].r
            else:
                r_b = np.array([body['orbital_radius'], 0])
                r_b1 = r_b - r_1
                vy = np.sqrt(self.G*m_1/norm(r_b1))
                self.bodies.append(Body(
                    body['name'], 
                    body['mass'],
                    body['orbital_radius'],
                    vy,
                    body['colour'],
                    body['display_radius']))
                               
                
    def TotalEnergy(self):
        """
        Computes the total energy (kinetic + gravitational potential) of the system.

        """
        KE_tot = 0.0
        GPE = 0.0
        for body in self.bodies:
            KE_tot += body.KineticEnergy()
            GPE += -1/2 * body.GPE(self.bodies, self.G)
        E_tot = GPE + KE_tot
        return E_tot
                
    def AnimateBodies(self, i):
        """
        Updates the positions of celestial bodies for animation.

        """
        self.update() #updates all positions, velocities and accelerations of bodies
        for c in range(len(self.patches)):
            self.patches[c].center = self.bodies[c].r #updates positions of the circles on the animated plot
        time = self.frame_number*self.dt
        self.time_text.set_text(f"Time: {round(time, 2)} years") #writes time on plot
        for body in self.bodies:
            body.period(self.bodies[0], time) #calls period function for each body
        return self.patches + [self.time_text]
        
    def display(self):
        """
        Sets up and displays the animation of celestial bodies.
        
        """
        fig = plt.figure()
        ax = plt.axes()
        
        self.patches = [] #empty list: it will contain position of all the circles representing celestial bodies
        
        #adds a point to the above list for each body
        for body in self.bodies: 
            point = plt.Circle(body.r, body.display_radius, color = body.colour, animated = True, label = body.name)
            self.patches.append(point)
            ax.add_patch(point)
        
        self.time_text = plt.text(0.4, 0.9, f"Time: {0} years", fontsize = 20, transform=ax.transAxes)
        
        ax.axis("scaled")
        ax.set_xlim(-6, 6) #can be changed to about -30, 30 to see outer planets
        ax.set_ylim(-6, 6) #can be changed to about -30, 30 to see outer planets
        ax.set_xlabel("X-Coordinate of Planet (au)")
        ax.set_ylabel("Y-Coordinate of Planet (au)")
        ax.legend(handles = self.patches)
        
        self.anim = FuncAnimation(fig, self.AnimateBodies, self.no_steps, repeat = False, interval = 10,
                                  blit = True) #animates bodies
        
        plt.show()
        
    def AnimateSatellites(self, i):
        """
        Updates the positions of celestial bodies and satellites for animation.

         """
        self.update() #updates positions of planets, sun AND satellite
        for c in range(len(self.patches)):
             self.patches[c].center = self.bodies[c].r #updates planet and sun patches to be plotted
        for d in range(len(self.dots)):
            self.dots[d].center =  self.satellites[d].r #updates satellite patches to be plotted
        time = self.frame_number*self.dt
        self.time_text.set_text(f"Time: {round(time, 2)} years")
        self.line.set_ydata(self.satellite_ypositions) #trajectory of satellite
        self.line.set_xdata(self.satellite_xpositions) #trajectory of satellite
        return self.patches + [self.time_text] + self.dots + [self.line,]
        
    def displaySatellite(self):
         """
         Sets up and displays the animation of satellites along with celestial bodies.
         
         """
         fig = plt.figure()
         ax = plt.axes()
         
         self.line, = ax.plot([],[], color = "black")
         
         self.patches = []
         
         for body in self.bodies:
             point = plt.Circle(body.r, body.display_radius, color = body.colour, animated = True, label = body.name)
             self.patches.append(point)
             ax.add_patch(point)
         
#The code that has been commented out below is the "Search Code". I used it to find a range of speeds to give my satellites before settling on one speed.
#I have left this code so it can be easily changed if I set a different goal for the satellite. 
         # satellite_speeds = np.linspace(0.87, 0.98, 10)
         # print(satellite_speeds)
         self.satellites = []
         # for v in satellite_speeds:    
         #     self.satellites.append(Satellite(self.bodies[3].r+(1/1000), (v + self.bodies[3].vy)))
         self.satellites.append(Satellite(self.bodies[3].r+(1/1000), (0.89 + self.bodies[3].vy)))
         self.dots = []
         
         for i in self.satellites: #there is only one item in this list but i left this loop so code is scalable
             dot = plt.Circle(i.r, 0.05, color = "black", animated = True, label = "Satellite")
             self.dots.append(dot)
             ax.add_patch(dot)
         
         self.time_text = plt.text(0.4, 0.9, f"Time: {0} years", fontsize = 20, transform=ax.transAxes)
         
         ax.axis("scaled")
         ax.set_xlim(-6, 6)
         ax.set_ylim(-6, 6)
         ax.set_xlabel("X-Coordinate of Planet (au)")
         ax.set_ylabel("Y-Coordinate of Planet (au)")
         ax.legend(handles = self.patches + self.dots)
         
         self.anim = FuncAnimation(fig, self.AnimateSatellites, self.no_steps, repeat = False, interval = 10,
                                   blit = True) #animates bodies and satellites
         
         plt.show()
        
    def WriteEnergies(self):
        """
        Updates the total energy over time and writes it to a file to be plotted.

        """
        def write_file(filename):
            with open(filename, "w") as file:
                for i in range(0,len(self.total_energy)):
                    file.write(f"{self.time[i]},{self.total_energy[i]}\n") #writes one point of time and energy per line in the file, separates them with a comma
        self.total_energy = []
        self.time = []
        for i in range(0,self.no_steps): #updates energy and time and adds them to a list for the number of iterations given in JSON file
            self.update()
            if i % 100 == 0: #only adds multiples of 100 to save time
                self.total_energy.append(self.TotalEnergy())
                self.time.append(i*self.dt)
        #identifies which method is being used and writes to the correct file        
        if self.method == "Beeman":       
            write_file("beeman.txt")
        elif self.method == "EulerCromer":
            write_file("eulercromer.txt")
        elif self.method == "Euler":
            write_file("euler.txt")

    
    def graph(self):
        """
        Displays a graph of total energy over time for the system.
        
        """
        def read_file(filename):
            with open(filename, "r") as file:
                lines = file.read().splitlines() #splits the file into individual lines
                for line in lines:
                    x, y = line.split(",") #dvides each line by the comma so distinct data can be accessed 
                    self.x_data.append(float(x)) #adds time to a list
                    self.y_data.append(float(y)) #adds energies
        self.x_data = []
        self.y_data = []
        #identifies which method is being used and reads the correct file
        if self.method == "Beeman":
            read_file("beeman.txt")
        elif self.method == "EulerCromer":
            read_file("eulercromer.txt")
        elif self.method == "Euler":
            read_file("euler.txt")
        line,= plt.plot(self.x_data[1:],self.y_data[1:]) #plots graph after first data point to avoid scaling issues
        plt.title("Total Energy of the System Vs Time")
        plt.xlabel("Time (earth years)")
        plt.ylabel("Total energy (earth units)")
        line.set_label(self.method)
        plt.legend()
        plt.show

class Beeman(Simulation):
    def __init__(self):
        super().__init__() 
        self.method = "Beeman" 
    def update(self):
        """
        Updates the acceleration, velocity, and position of bodies and satellites 
        using Beeman's integration method.

        The method ensures accurate position updates by considering the previous 
        acceleration values for better energy conservation.
        
        Notes
        -----
        This function uses the Euler-Cromer method for the first timestep as it does 
        not require an acceleration at the previous timestep.
        For the second time step, there is no value for previous acceleration so it 
        use the current acceleration of the body.
            
        """
        if self.frame_number == 0:
            for body in self.bodies: #updates all accelerations of the planets before updating positions and velocities 
                body.a = body.accelerate(self.bodies, self.G)
            for body in self.bodies:
                body.v = body.v + body.a*self.dt
                body.r = body.r + body.v*self.dt
        else:
            for body in self.bodies:
                if body.a_prev is None:
                    body.a_prev = body.a #uses current accel if prev accel has no value
                body.r = body.r + body.v*self.dt + (1/6)*(4*body.a - body.a_prev)*(self.dt)**2
                
            for body in self.bodies:
                a_prev2 = body.a_prev #updates names in preparation of next time step
                body.a_prev = body.a 
                body.a = body.accelerate(self.bodies, self.G)
                body.v = body.v + (1/6)*(2*body.a + 5*body.a_prev - a_prev2)*self.dt
        
        if self.satellites:
            e_and_m = self.bodies[3:5] #earth and mars bodies 
            if self.frame_number == 0: #euler-cromer method for first timestep
                for s in self.satellites:
                    s.a = s.accelerate(self.bodies, self.G)
                for s in self.satellites:
                    s.v = s.v + s.a*self.dt
                    s.r = s.r + s.v*self.dt
            else: 
                for s in self.satellites: #beeman method
                    if s.a_prev is None:
                        s.a_prev = s.a
                    s.r = s.r + s.v*self.dt + (1/6)*(4*s.a - s.a_prev)*(self.dt)**2
                    s.closestDistance(e_and_m, self.frame_number*self.dt)
                    
                    
                for s in self.satellites:
                    a_prev2 = s.a_prev
                    s.a_prev = s.a 
                    s.a = s.accelerate(self.bodies, self.G)
                    s.v = s.v + (1/6)*(2*s.a + 5*s.a_prev - a_prev2)*self.dt
            
        
            self.satellite_xpositions.append(self.satellites[0].r[0])
            self.satellite_ypositions.append(self.satellites[0].r[1])
        self.frame_number += 1 #update frame number which can be used to calculate time in other parts of program
        
class EulerCromer(Simulation):
    def __init__(self):
        super().__init__()
        self.method = "EulerCromer"
    def update(self):
        """
        Updates the acceleration, velocity, and position of bodies using the Euler-Cromer method.
            
        """
        for body in self.bodies: #updates all accelerations first before updating velocities and positions
            body.a = body.accelerate(self.bodies, self.G)
            
        for body in self.bodies:
            body.v = body.v + body.a*self.dt
            body.r = body.r + body.v*self.dt
        
        
class Euler(Simulation):
    def __init__(self):
        super().__init__()
        self.method = "Euler"
    def update(self):
        """
        Updates the acceleration, velocity, and position of bodies using the Direct Euler integration method.
            
        """
        for body in self.bodies:
            body.a = body.accelerate(self.bodies, self.G)
            
        for body in self.bodies:
            body.r = body.r + body.v*self.dt
            body.v = body.v + body.a*self.dt
    
    
def main():
    """
    Initialises and runs simulations using different numerical integration methods.

    This function:
    - Creates instances of the `Beeman`, `EulerCromer`, and `Euler` simulation classes.
    - Displays the simulation without a satellite using Beeman integration method
    - Displays the satellite animation using the Beeman integration method.
    - Displays energy graphs for comparison.
    """
    sim1 = Beeman()
    sim2 = EulerCromer()
    sim3 = Euler()
    #Experiment 1: Orbits
    sim1.display()
    
    #Experiment 3: Satellite
    #sim1.displaySatellite()
    
    #Experiment 2: Energies
    #sim1.WriteEnergies()
    #sim1.graph()
    #sim2.WriteEnergies()
    #sim2.graph()
    #sim3.WriteEnergies()
    #sim3.graph()
main()