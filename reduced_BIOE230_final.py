from random import random
import math
import matplotlib.pyplot as plt
from matplotlib import colors

# ABP model parameters: ballistic velocity, time step, translational diffusion constant
velocity_magnitude = 8
dt = 1
rotational_diffusion = 0.15
translational_diffusion = 5
#vel =5.0; dt = 1; Drot = 0.1; Dtrans = 0.1; #original parameters that work well in case the new ones get screwed up

#Initialize scaling stuff (how long the sim will take)
max_steps = 10000             #will remove by having termination functionality
num_walkers = 15
x_scaling = 2000            #width of area walkers could be on

#Nutrient info. Pretty sure you can change the "center" value
nutrient_level = 1500       # point where nutrient is located
nutrient_exponent = 1.5       # Concentration is 10^exponent (force strength)
start_exponent = 0.5          # Concentration is 10^exponent (force strength)
# nutrient_exponent, start_exponent = 5, 3 # Concentration is 10^exponent (force strength with -6)


class Walker():
    def __init__(self):
        self.starting_position = [x_scaling*random(), 0]
        self.current_position = self.starting_position
        self.terminated = False
        self.positions = [self.current_position]
        self.velocity_direction = 2*math.pi*random()

    def walk(self):
        '''
        executes one step of walk and updates position
        '''

        #calculate force and random vectors and update position
        force_vector = self.gradient_force()
        random_vector = self.random_force()
        total_vector = [force_vector[0] + random_vector[0], force_vector[1] + random_vector[1]]
        total_vector_length = math.sqrt(total_vector[0]**2 + total_vector[1]**2)
        total_vector = [dt*total_vector[0]/total_vector_length, dt*total_vector[1]/total_vector_length]
        self.current_position = [self.current_position[0] + total_vector[0], self.current_position[1] + total_vector[1]]
        self.positions.append(self.current_position)

        self.terminate()
        self.branch()

    def gradient_force(self):
        '''
        Calculates contribution of contribution gradient on movement
        '''        
        distance_to_nutrient = abs(self.current_position[1] - nutrient_level)
        noramlized_distance_to_nutrient = 1 - (distance_to_nutrient/nutrient_level)
        exponent = (noramlized_distance_to_nutrient * (nutrient_exponent - start_exponent)) + start_exponent
        y_force = 10 ** exponent

        return [0, y_force]

    def random_force(self):
        '''
        Calculates contrbution of randomness on movement
        '''
        #calculate translational randomness
        dx = translational_diffusion * dt * 2*(random() - 0.5)
        dy = translational_diffusion * dt * 2*(random() - 0.5)

        #calculate new velocity
        velocity = [velocity_magnitude*dt*math.cos(self.velocity_direction), velocity_magnitude*dt*math.sin(self.velocity_direction)]

        #update valocity direction with some randomness
        self.velocity_direction += rotational_diffusion * dt * (2*math.pi*(random()-0.5))

        return [velocity[0] + dx, velocity[1] + dy]

    def terminate(self):
        '''
        Used to determine if the walker should be terminated or is finished and executes the termination
        '''
        if self.current_position[1] > nutrient_level:
            self.terminated = True

    def branch(self):
        '''
        Used to determine if the walker should branch and executes the branch
        '''
        pass

def ABP_step_biased_linear(plot = True):
    '''
    Runs walker simulation and plots results
    '''
    #generate initial list of walkers
    walkers = [Walker() for i in range(num_walkers)]
    finished_walkers = []

    #iterate through steps
    steps = 0
    while len(walkers) > 0 and steps < max_steps:
        for walker in walkers:
            walker.walk()
            if walker.terminated:
                walkers.remove(walker)
                finished_walkers.append(walker)
        steps += 1

    if plot:
        fig, ax = plt.subplots(1, 1, figsize = (10, 10))
        all_walkers = walkers + finished_walkers
        for walker in all_walkers:
            #print(walker.positions)
            x = []
            y = []
            for position in walker.positions:
                x.append(position[0])
                y.append(position[1])
            plt.plot(x, y, '-')
            plt.plot(x[-1], y[-1], 'ro', markersize = 8)
            plt.plot(x[0], y[0], 'bo', markersize = 8)
        ax.set_title("Biased walk with chemoattractant".format(), x = 0.5, y = 0.87)
        ax.set_xlim(-1000, 3000)
        ax.set_ylim(-1000, 3000)
        ax.set_xlabel("poisiton in μm")
        ax.set_ylabel("poisiton in μm")

        #First set color map
        mycolor = [[256, 256, 256], [256, 255, 254], [256, 253, 250], [256, 250, 240], [255, 236, 209], [255, 218, 185], [251, 196, 171], [248, 173, 157], [244, 151, 142], [240, 128, 128]] #from coolors：）
        for i in mycolor:
            for j in range(len(i)):
                i[j] *= (1/256)
        cmap_color = colors.LinearSegmentedColormap.from_list('my_list', mycolor)
        m = 4000
        conc_matrix = [0]*m    #np.zeros((m, m))
        for i in range(m):
            distance_to_nutrient = abs(i - 1000 - nutrient_level)
            noramlized_distance_to_nutrient = 1 - (distance_to_nutrient/nutrient_level)
            exponent = (noramlized_distance_to_nutrient * (nutrient_exponent - start_exponent)) + start_exponent
            conc_matrix[i] = [exponent + 6]*m
        ax.imshow(conc_matrix, cmap=cmap_color, interpolation='nearest', extent = [-1000, 3000, -1000, 3000], origin = 'lower')

        plt.show()


if __name__ == '__main__':
    '''
    run code if you run this script instead of importing it
    '''
    ABP_step_biased_linear()