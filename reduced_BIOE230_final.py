from random import random
import math
import matplotlib.pyplot as plt
from matplotlib import colors

# ABP model parameters: ballistic velocity, time step, translational diffusion constant
                

#Initialize scaling stuff (how long the sim will take)
num_walkers = 5
x_scaling = 2000            #wiself.dth of area walkers could be on

#Nutrient info. Pretty sure you can change the "center" value
nutrient_level = 1500       # point where nutrient is located
# self.nutrient_exponent, self.start_exponent = 5, 3 # Concentration is 10^exponent (force strength with -6)



'''
steps move to walker so new walkers start at step 0
    what can be added to walker:
        max_steps
        magnitudes
        self.dt
        diffusionsexponents
        self.branch_probability
        con/div stuff

add break command that plots whatever is currently computed

update print outs to say how far along we are

run time, maybe a log saying what the variables were as well

save plots to folder

make variables seem like they mean something

convergence_divergence
    wait longer to converge?
    maybe x^3 function with 0.5 triple root

'''


class Walker():
    def __init__(self, starting_position = None):
        '''
        define constants and initialized internal variables
        '''
        #define constants
        #max steps and step size
        self.max_steps = 3000
        self.dt = 1

        #magnitudes
        self.momentum_magnitude = 15                  
        self.random_magnitude = 5
        self.convergence_divergence_magnitude = 25

        #behavioral constants
        self.rotational_diffusion = 0.05
        self.nutrient_exponent = 1.5       # Concentration is 10^exponent (force strength)
        self.start_exponent = 1        # Concentration is 10^exponent (force strength)
        self.branch_probability = 0.0005
        self.convergence_divergence_distance = 100

        #initialize variables
        if starting_position == None:
            self.starting_position = [x_scaling*random(), 0]
        else:
            self.starting_position = starting_position
        self.current_position = self.starting_position
        self.finished = False
        self.positions = [self.current_position]
        self.velocity_direction = 2*math.pi*random()
        self.branched = False
        self.convergence_divergence_vector = [0, 0]
        self.terminated = False
        self.step = 0

    def walk(self, walkers):
        '''
        executes one step of walk and updates position
        '''

        #calculate force and random vectors and update position
        force_vector = self.gradient_force()
        random_vector = self.random_force()
        total_vector = [force_vector[0] + random_vector[0] + self.convergence_divergence_vector[0], force_vector[1] + random_vector[1] + self.convergence_divergence_vector[1]]
        total_vector_length = math.sqrt(total_vector[0]**2 + total_vector[1]**2)
        total_vector = [self.dt*total_vector[0]/total_vector_length, self.dt*total_vector[1]/total_vector_length]
        self.current_position = [self.current_position[0] + total_vector[0], self.current_position[1] + total_vector[1]]
        self.positions.append(self.current_position)

        self.finish()
        if self.finished:
            return
        self.convergence_divergence_vector = [0, 0]
        for comparison_walker in walkers:
            self.check_terminate_and_update_convergence_divergence_force(comparison_walker)
            if self.terminated:
                break

        self.step += 1

    def gradient_force(self):
        '''
        Calculates contribution of contribution gradient on movement
        '''        
        distance_to_nutrient = abs(self.current_position[1] - nutrient_level)
        noramlized_distance_to_nutrient = 1 - (distance_to_nutrient/nutrient_level)
        exponent = (noramlized_distance_to_nutrient * (self.nutrient_exponent - self.start_exponent)) + self.start_exponent
        y_force = 10 ** exponent

        return [0, y_force]

    def random_force(self):
        '''
        Calculates contrbution of randomness on movement
        '''
        #calculate translational randomness
        dx = self.random_magnitude * 2 *(random() - 0.5)
        dy = self.random_magnitude * 2 *(random() - 0.5)

        #calculate new velocity
        velocity = [self.momentum_magnitude*math.cos(self.velocity_direction), self.momentum_magnitude*math.sin(self.velocity_direction)]

        #update valocity direction with some randomness
        self.velocity_direction += self.rotational_diffusion * (2*math.pi*(random()-0.5))

        return [velocity[0] + dx, velocity[1] + dy]

    def check_collision(self, previous_comparison_position, comparison_position):
        '''
        checks to see if there is a collision between two segments
        '''
        '''
        next_segment_flag_1 = False
        do_break = False
        previous_position = self.positions[-2]
        lower_y = round(min(previous_position[1], self.current_position[1]), 5)
        upper_y = round(max(previous_position[1], self.current_position[1]), 5)
        lower_comparison_y = round(min(previous_comparison_position[1], comparison_position[1]), 5)
        upper_comparison_y = round(max(previous_comparison_position[1], comparison_position[1]), 5)

        if (lower_y < round(comparison_position[1], 5) < upper_y):
            next_segment_flag_1 = True
            #run check on previous segement
            y_solution = self.get_interception(comparison_position, previous_comparison_position)
            if (max(lower_y, lower_comparison_y) < y_solution < min(upper_y, upper_comparison_y)):
                self.terminated = True
                do_break = True
        elif next_segment_flag:
            #run check on previous segement
            y_solution = self.get_interception(comparison_position, previous_comparison_position)
            if (max(lower_y, lower_comparison_y) < y_solution < min(upper_y, upper_comparison_y)):
                self.terminated = True
                do_break = True
        elif lower_comparison_y < lower_y and upper_comparison_y > upper_y:
            #run check
            y_solution = self.get_interception(comparison_position, previous_comparison_position)
            if (lower_y < y_solution < upper_y):
                self.terminated = True
                do_break = True

        return next_segment_flag_1, do_break
        '''
        do_break = False
        y_solution = self.get_interception(comparison_position, previous_comparison_position)
        lower_y = round(max(min(self.positions[-2][1], self.current_position[1]), min(previous_comparison_position[1], comparison_position[1])), 5)
        upper_y = round(min(max(self.positions[-2][1], self.current_position[1]), max(previous_comparison_position[1], comparison_position[1])), 5)
        if lower_y < y_solution < upper_y:
            self.terminated = True
            do_break = True

        return do_break


    def check_terminate_and_update_convergence_divergence_force(self, comparison_walker):
        '''
        loops through other walkers and avoids them for a bit and then goes towards them as time goes on
        '''
        normalized_time = (self.step/self.max_steps) - 0.5
        previous_comparison_position = comparison_walker.starting_position
        next_segment_flag = False

        skip_first = True
        for comparison_position in comparison_walker.positions:
            #iterate through positions on the comparison walker

            if skip_first:
                skip_first = False
                continue

            distance_to_position = round(math.sqrt((comparison_position[0]-self.current_position[0])**2 + (comparison_position[1]-self.current_position[1])**2), 5)
            if distance_to_position == 0:
                if comparison_walker == self:
                    continue
                self.terminated = True
                break

            if distance_to_position < self.dt*2:
                do_break = self.check_collision(previous_comparison_position, comparison_position)
                if do_break:
                    break

            previous_comparison_position = comparison_position

            if distance_to_position < self.convergence_divergence_distance:
                normalized_vector_magnitude = self.convergence_divergence_magnitude * normalized_time / distance_to_position**2
                self.convergence_divergence_vector[0] += normalized_vector_magnitude * (comparison_position[0]-self.current_position[0])
                self.convergence_divergence_vector[1] += normalized_vector_magnitude * (comparison_position[1]-self.current_position[1])

    def get_interception(self, comparison_position, previous_comparison_position):
        previous_position = self.positions[-2]

        current_slope = (previous_position[0] - self.current_position[0])/(previous_position[1] - self.current_position[1])
        comparison_slope = (previous_comparison_position[0] - comparison_position[0])/(previous_comparison_position[1] - comparison_position[1])
        numerator = previous_position[1]*current_slope + previous_comparison_position[0] - previous_position[0] - previous_comparison_position[1]*comparison_slope
        denominator = current_slope - comparison_slope
        return round(numerator/denominator, 5)


    def finish(self):
        '''
        Used to determine if the walker is finished
        '''
        if self.current_position[1] > nutrient_level or self.step > self.max_steps:
            self.finished = True

    def branch(self):
        '''
        Used to determine if the walker should branch and executes the branch
        '''
        if random() < self.branch_probability:
            return True
        return False

def ABP_step_biased_linear(plot = True):
    '''
    Runs walker simulation and plots results
    '''
    '''
    #generate initial list of walkers
    walkers = [Walker() for i in range(num_walkers)]
    finished_walkers = []

    #iterate through steps
    steps = 0
    while len(walkers) > 0 and steps < max_steps:
        for walker in walkers:
            walker.walk(walkers+finished_walkers)
            if walker.finished:
                walkers.remove(walker)
                finished_walkers.append(walker)
                continue
            for comparison_walker in walkers+finished_walkers:
                #if comparison_walker == walker:
                    #continue
                walker.check_terminate_and_update_convergence_divergence_force(steps, comparison_walker)
                if walker.terminated:
                    walkers.remove(walker)
                    finished_walkers.append(walker)
                    break
            if walker.terminated:
                continue
            if walker.branch():
                walkers.append(Walker(walker.current_position))
        steps += 1
    print(steps)
    '''
    walkers = [Walker() for i in range(num_walkers)]
    finished_walkers = []

    while len(walkers) > 0:
        for walker in walkers:
            walker.walk(walkers+finished_walkers)
            if walker.finished or walker.terminated:
                walkers.remove(walker)
                finished_walkers.append(walker)
                continue
            if walker.branch():
                walkers.append(Walker(walker.current_position))




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
            plt.plot(x[-1], y[-1], 'ro', markersize = 4)
            plt.plot(x[0], y[0], 'bo', markersize = 4)
        ax.set_title("Biased walk with chemoattractant".format(), x = 0.5, y = 0.87)
        ax.set_xlim(-3000, 5000)
        ax.set_ylim(-3000, 3000)
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
            exponent = noramlized_distance_to_nutrient * 2 #(noramlized_distance_to_nutrient * (self.nutrient_exponent - self.start_exponent)) + self.start_exponent
            conc_matrix[i] = [exponent + 6]*m
        ax.imshow(conc_matrix, cmap=cmap_color, interpolation='nearest', extent = [-1000, 3000, -1000, 3000], origin = 'lower')

        plt.show()


if __name__ == '__main__':
    '''
    run code if you run this script instead of importing it
    '''
    ABP_step_biased_linear()