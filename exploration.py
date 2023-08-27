from random import random
import math
import matplotlib.pyplot as plt
from matplotlib import colors
import msvcrt
import time
from datetime import datetime
import os

# Simulation Constants (Walker settings in Walker.__init__ function)
num_walkers = 5             # Number of Walkers
x_scaling = 2000            # width of area walkers could be on
nutrient_level = 1500       # point where nutrient is located
max_loops = 10000           # stop the simulation from running away


class Walker():
    def __init__(self, starting_position = None):
        '''
        define constants and initialized internal variables
        '''
        #define constants
        #max steps and step size
        self.max_steps = 3000                                       #maximum steps each walker can take
        self.dt = 1                                                 #step size

        #magnitudes
        self.momentum_magnitude = 15                                #relative magnitude of the momemtum vector
        self.random_magnitude = 5                                   #relative magnitude of the random vector
        self.convergence_divergence_magnitude = 25                  #relative magnitude of the convergence/diveregence vector
        self.nutrient_exponent = 1.5                                #force vector magnitude is 10^[nutrient_exponent] at the nutrient
        self.start_exponent = 1                                     #force vector magnitude is 10^[start_exponent] at the start
                                                                    #force vetor expoenet changes linearly as the walker approaches the nutrient
        #behavioral constants
        self.rotational_diffusion = 0.05                            #how much the momentum can change per step
        self.branch_probability = 0.0005                            #probability that a walker with branch
        self.convergence_divergence_distance = 100                  #radius the walker will look to determine the convergence/divergence force

        #initialize variables (Don't Touch!)
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
        force_vector = self.gradient_force()                    #calculate force vector
        random_vector = self.random_force()                     #calculate random vector
        momentum_vector = self.momentum_force()                      #calculate momentum vector

        #Sum up all vectors (force, random, momentum, convergence/divergence)
        total_vector = [force_vector[0] + 
                        random_vector[0] + 
                        momentum_vector[0] + 
                        self.convergence_divergence_vector[0], 
                        force_vector[1] + 
                        random_vector[1] + 
                        momentum_vector[1] + 
                        self.convergence_divergence_vector[1]]

        #calculate vecotr length to normailze
        total_vector_length = math.sqrt(total_vector[0]**2 + total_vector[1]**2)

        #normalize vector length to make total vector of length dt
        total_vector = [self.dt*total_vector[0]/total_vector_length, 
                        self.dt*total_vector[1]/total_vector_length]

        #update current position
        self.current_position = [self.current_position[0] + total_vector[0], 
                                 self.current_position[1] + total_vector[1]]
        self.positions.append(self.current_position)

        #Check if walker reached the nutrient or max_steps
        self.finish()
        if self.finished or self.terminated:
            return

        #check if the walker ran into another walker. Terminate if it did
        #otherwise, update convergence/divergence vector for next step
        #Loops combined for efficiency
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
        dx = self.random_magnitude * 2 *(random() - 0.5)
        dy = self.random_magnitude * 2 *(random() - 0.5)

        return [dx, dy]

    def momentum_force(self):
        '''
        Calculates contribution of momentum on movement
        '''
        #update velocity direction with some randomness
        self.velocity_direction += self.rotational_diffusion * (2*math.pi*(random()-0.5))

        #calculate new velocity
        velocity = [self.momentum_magnitude*math.cos(self.velocity_direction), self.momentum_magnitude*math.sin(self.velocity_direction)]

        return velocity


    def check_collision(self, previous_comparison_position, comparison_position, distance_to_position):
        '''
        checks to see if there is a collision between two segments
        '''
        do_break = False

        #Only check collisions if the comparison position is within 2 segment lenghts
        if distance_to_position < self.dt*2:
            y_solution = self.get_interception(comparison_position, previous_comparison_position)
            lower_y = round(max(min(self.positions[-2][1], self.current_position[1]), min(previous_comparison_position[1], comparison_position[1])), 5)
            upper_y = round(min(max(self.positions[-2][1], self.current_position[1]), max(previous_comparison_position[1], comparison_position[1])), 5)
            if lower_y < y_solution < upper_y:
                self.positions[-1] = comparison_position
                self.terminated = True
                do_break = True

        return do_break


    def check_terminate_and_update_convergence_divergence_force(self, comparison_walker):
        '''
        loops through other walker's positions
        checks for collisions for termination
        updates convergence/diverence if no collisions
        '''
        previous_comparison_position = comparison_walker.starting_position      #initialize
        normalized_time = ((self.step/self.max_steps) - 0.5)**3                 #used to determine convergence v. divergence and part of the magnitude of the force (cubic function)

        skip_first = True
        for comparison_position in comparison_walker.positions:
            #iterate through positions on the comparison walker
            if skip_first:
                #need first segment not first point
                skip_first = False
                continue

            #calculate distance to comparison point
            distance_to_position = round(math.sqrt((comparison_position[0]-self.current_position[0])**2 + (comparison_position[1]-self.current_position[1])**2), 5)
            
            if distance_to_position == 0:
                #can only really happen if they are the same point. Ignoring this case
                continue

            #check for collision between new segment and comparison segment
            do_break = self.check_collision(previous_comparison_position, comparison_position, distance_to_position)
            if do_break:
                break

            #update comparison position for next loop
            previous_comparison_position = comparison_position

            #calculate convergence/diverence force for next step if comparision position 
            if distance_to_position < self.convergence_divergence_distance:
                #normalize force vector for controlled magnitude
                normalized_vector_magnitude = self.convergence_divergence_magnitude * normalized_time / distance_to_position**2

                #update convergence divergence vector
                self.convergence_divergence_vector[0] += normalized_vector_magnitude * (comparison_position[0]-self.current_position[0])
                self.convergence_divergence_vector[1] += normalized_vector_magnitude * (comparison_position[1]-self.current_position[1])

    def get_interception(self, comparison_position, previous_comparison_position):
        '''
        Find y value of the intersection between the extrapolated segements
        '''
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
        if self.current_position[1] > nutrient_level:
            self.finished = True
        if self.step > self.max_steps:
            self.terminated = True

    def branch(self):
        '''
        Used to determine if the walker should branch and executes the branch
        '''
        if random() < self.branch_probability:
            return True
        return False

def exploitation_algorithm():
    '''
    Runs walker simulation and plots results
    '''
    #Instructions
    print('Hit Enter if you want to stop the simulation early. It still may take a second to stop')

    start_time = time.time()                                        #start run timer

    #Initialize walker lists
    walkers = [Walker() for i in range(num_walkers)]
    finished_walkers = []

    loop_count = 0
    while len(walkers) > 0 and loop_count < max_loops:
        #loop simulation while there are walkers not done yet
        for walker in walkers:
            #loop through walkers so they step together
            walker.walk(walkers+finished_walkers)                   #execute step
            if walker.finished or walker.terminated:                #move wakler to finished list of it finished or terminated
                walkers.remove(walker)
                finished_walkers.append(walker)
                continue
            if walker.branch():                                     #create new walker if branched
                walkers.append(Walker(walker.current_position))

        #protect run-away loop and print progress
        loop_count += 1
        if loop_count%100 == 0:
            print('Loop count: ', loop_count)
            print('Walkers remaining: ', len(walkers))

        #end simulation early if enter is pressed
        if msvcrt.kbhit():
            if msvcrt.getwche() == '\r':
                break

    end_time = time.time()
    print(f'Simulation time: {round(end_time - start_time, 0)} seconds')

    #log results
    with open('log.csv', 'a') as log_file:
        '''
        save simulation parameters to log
        date,time,number of walkers,nutrient level,number of total walkers,number of finished walkers,
        max_steps,dt,momentum magnitude,random magnitude,convergence/divergence magnitude,
        nutrient exponent,start exponent,rotational diffusion,branch probability,
        convergence/divergence distance,loop count,simulation time
        '''
        now = datetime.now()                #get date for saves and logs
        log_file.write(f'{now.strftime("%m/%d/%Y")},{now.strftime("%H:%M:%S")},')    #enter date and time
        log_file.write(f'{num_walkers},')
        log_file.write(f'{nutrient_level},')
        log_file.write(f'{len(walkers+finished_walkers)},')
        finshed_walker_count = 0
        for walker in finished_walkers:
            if walker.finished:
                finshed_walker_count += 1
        log_file.write(f'{finshed_walker_count},')
        all_walkers = walkers+finished_walkers
        log_file.write(f'{all_walkers[0].max_steps},')
        log_file.write(f'{all_walkers[0].dt},')
        log_file.write(f'{all_walkers[0].momentum_magnitude},')
        log_file.write(f'{all_walkers[0].random_magnitude},')
        log_file.write(f'{all_walkers[0].convergence_divergence_magnitude},')
        log_file.write(f'{all_walkers[0].nutrient_exponent},')
        log_file.write(f'{all_walkers[0].start_exponent},')
        log_file.write(f'{all_walkers[0].rotational_diffusion},')
        log_file.write(f'{all_walkers[0].branch_probability},')
        log_file.write(f'{all_walkers[0].convergence_divergence_distance},')
        log_file.write(f'{loop_count},')
        log_file.write(f'{end_time - start_time},')
        log_file.write(f'{input("Enter a note for the log file: ")}\n')


    return walkers + finished_walkers







def save_data(all_walkers, plot=None):
    now = datetime.now()                #get date for saves and logs
    #save results
    #make new directory
    os.mkdir(f'raw_data\\{now.strftime("%m_%d_%Y")}_{now.strftime("%H-%M-%S")}')

    if plot != None:
        #save plot
        plot.savefig(f'raw_data\\{now.strftime("%m_%d_%Y")}_{now.strftime("%H-%M-%S")}\\plot.png')

    #generate and save walker data
    with open(f'raw_data\\{now.strftime("%m_%d_%Y")}_{now.strftime("%H-%M-%S")}\\raw_data.csv', 'a') as walker_file:
        walker_file.write('walker_id,position_id,x,y\n')
        walker_id = 0
        for walker in all_walkers:
            position_id = 0
            for position in walker.positions:
                walker_file.write(f'{walker_id},{position_id},{position[0]},{position[1]}\n')
                position_id += 1
            walker_id += 1
    
def plot_data(all_walkers):
    #plot results
    fig, ax = plt.subplots(1, 1, figsize = (10, 10))
    for walker in all_walkers:
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
    m = 8000
    conc_matrix = [0]*m    #np.zeros((m, m))
    for i in range(m):
        distance_to_nutrient = abs(i - 3000 - nutrient_level)
        noramlized_distance_to_nutrient = 1 - (distance_to_nutrient/nutrient_level)
        exponent = noramlized_distance_to_nutrient * 2 #(noramlized_distance_to_nutrient * (self.nutrient_exponent - self.start_exponent)) + self.start_exponent
        conc_matrix[i] = [exponent + 6]*m
    ax.imshow(conc_matrix, cmap=cmap_color, interpolation='nearest', extent = [-3000, 5000, -3000, 5000], origin = 'lower')

    return plt

    
if __name__ == '__main__':
    '''
    run code if you run this script instead of importing it
    '''
    all_walkers = exploitation_algorithm()
    save_data(all_walkers)
    plot_data(all_walkers)