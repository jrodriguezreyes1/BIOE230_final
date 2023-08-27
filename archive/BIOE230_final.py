import numpy as np
import matplotlib.pyplot as plt
from numpy.random import random as rand
from scipy import stats
import math
from matplotlib import colors
import os

'''
This will be the conglomerate of the ipynb file
'''

'''
3 main points
branching
recombining
divergence to convergence with time
'''

'''
I have a vision of how we can make this simulation a lot better.  
I feel like every singular line is straightforward, albeit computationally excessive overall.  
I'm not afraid to use the cluster.  
I know it won't be super easy to implement, but I can't seem to even initiate a damn thing.  
Let me know what you think of it.  
I totally value your opinion if you see somewhere it can be done more efficiently or just better.  
Also lmk if you think of a task that is a better use of my time than fighting with np.append.
Do the following right after EVERY STEP for EVERY WALKER (already a nested for loop)

temp_rand=rand(1)
If temp_rand > 0.95
	np.append(new_nodes, [x,y])

For n=len(master_coords)
#Calculate Euclidean distance between current step’s coordinates and n coordinates
	If distance <= 10
		end walk
	Elseif 
		x += c * (100 - num_step) * (1/(x - master_coords[1,n]))
	np.append(master_coords,[x,y])

In English:
Give it a 1/XX chance (exact chance will be determined later) to create a new walker. 
If it passed this test, save coordinates into new_nodes and proceed.
First random walker is initiated.  It takes one step determined by the code already made.   
Calculate distance to all past steps (master_coords) one at a time.  
If it hits itself (distance less than 10 or whatever), end the walk, ie, proceed to the next j.  
If this criterion is not met, run this calculation: 
	x += c * (100 - num_step) * (x - all_coords[1]).  
	“c”= fitting parameter to make it look right, 
	“100-num_steps” to make it repulsive/exploratory for a young walker and attractive/node seeking for an older walker, 
	“1/(x-master_coords[1,n])” to get direction and intensity of force (when x and the x component of the test coordinates are close, the force, either attractive or repulsive is strong, when this difference is large, the force is very low.” The final coordinates are added to the matrix master_coords.

The next walker does the same thing.  
More walkers= higher chance of hitting another coordinate and being cut off.  This is part of the plan.
After the 10 original walkers are done, then we move onto new_nodes.  
Treat each of these coordinates as we did the original 10.  
Some of these will be terminated, some will reach the end point, some will create more new nodes.  
Continue until there are no more new nodes.
'''

#  Objective 2.1: Model random walkers trajectories that move in response to chemoattractants scattered along a straight line# ABP model parameters: ballistic velocity, time step, rotational diffusion constant, translational diffusion constant

# Removed x term from distance calculation so only y vector of chemoattractant is considered




# Initialize parameters of ABP model and locations first

# ABP model parameters: ballistic velocity, time step, rotational diffusion constant, translational diffusion constant
vel =8.0; dt = 1; Drot = 0.03; Dtrans = 3;
#vel =5.0; dt = 1; Drot = 0.1; Dtrans = 0.1; #original parameters that work well in case the new ones get screwed up

# initialize arrays that store x,y and theta values, 
# as well as initial particle position and angle
num_steps = 800;
num_walks = 15;
x_scaling = 2000;

nutrient_center = [0, 1500] # point where nutrient is located
nutrient_exponent, start_exponent = 5, 3 # Concentration is 10^exponent for each respective one
# nutrient_exponent, start_exponent = 5, 3 # Concentration is 10^exponent for each respective one
proximity = 0 # distance from start to center, will be calculated later


# Calculates the distance from the random walkers current position (a) to the focal point of the nutrient (b)
def distance_to_nutrient(a,b):
    return math.sqrt(( (a[1] - b[1]) ** 2))


# Calculates concentration at current position
def calculate_conc(pos):
    dist = distance_to_nutrient(pos, nutrient_center) #
    initial_proximity = distance_to_nutrient([0,0],nutrient_center); #initialize proximity to nutrient at start
    exponent = (1 - dist / initial_proximity) * (nutrient_exponent - start_exponent) + start_exponent
    return 10 ** exponent


# Calculates contribution of contribution gradient on movement
def gradient_force(x,y,xnutr,ynutr):
    proximity = distance_to_nutrient([x,y],nutrient_center)
    
    concentration = calculate_conc([x,y])
    
    x_force = 0*(xnutr - x)*concentration*10**-6; # as concentration increases, make walker move towards 
    y_force = (ynutr - y)*concentration*10**-6; # as concentration increases, make walker move towards 
    
    x_combined = x+x_force;
    y_combined = y+y_force;
    return x_combined, y_combined, x_force, y_force, proximity, concentration

    # Define biased random walker with chemoattractant where steps are taken using ABP and chemoattractant is set across
# a line at the top


def ABP_step_biased_linear():

    for j  in range(num_walks):
        xvec=np.zeros(0); yvec=np.zeros(0); thetavec = np.zeros(0);
        x=x_scaling*(rand()); y = 0.0; theta = (2*np.pi)*rand();
        
        proximity = distance_to_nutrient([x,y],nutrient_center); #initialize proximity to nutrient at start
#         print('Initial proximity =',proximity)
        concentration = calculate_conc([x,y]) # initialize concentration of nutrient at start
#         print('Initial concentration =',concentration,'\n')
        
        # inner for loop is each step for a given walker/trajectory        
        for i in range(num_steps):
                # calculate diffusive/random steps. For the x- and y-,we generate 
                #a random number between -1 & 1 to ensure that the walker can step in both directions(up/down and left/right).
                
                # calculate change in position
                dx = np.sqrt(2*Dtrans*dt)*2*(rand(1)-0.5); 
                dy= np.sqrt(2*Dtrans*dt)*2*(rand(1)-0.5); 
                dtheta = np.sqrt(2*Drot*dt)*(2*np.pi)*(rand(1)-0.5);

                # update coordinates (including ballistic step)
                x += vel*dt*np.cos(theta) + dx;
                y += vel*dt*np.sin(theta) + dy;

                # store successive positions in arrays
                xvec = np.append(xvec,x); yvec = np.append(yvec,y); 
                # update the angle and store in array
                theta += dtheta;
                thetavec = np.append(thetavec, theta);
                

                # store successive positions in arrays
                pos_comb = gradient_force(x,y,nutrient_center[0],nutrient_center[1])
                x = pos_comb[0];
                y = pos_comb[1];
                                
        plt.plot(xvec,yvec, '-');
        plt.plot(xvec[-1],yvec[-1],'ro', markersize = 8)
        plt.plot(xvec[0],yvec[0],'bo', markersize = 8)
#         plt.plot(nutrient_center[0],nutrient_center[1],'g*', markersize = 16)

ABP_step_biased_linear()


# Make plot of chemoattractant gradient with biased random walk that responds to chemoattractant
plot_lim_x = [-1000,3000];
plot_lim_y = [-1000,3000];


#Below are all for plotting purposes
fig, ax = plt.subplots(1, 1, figsize = (10, 10))

#First set color map
mycolor = [[256, 256, 256], [256, 255, 254], [256, 253, 250], [256, 250, 240], [255, 236, 209], [255, 218, 185], [251, 196, 171], [248, 173, 157], [244, 151, 142], [240, 128, 128]] #from coolors：）
for i in mycolor:
    for j in range(len(i)):
        i[j] *= (1/256)
cmap_color = colors.LinearSegmentedColormap.from_list('my_list', mycolor) #Linearly segment these colors to create a continuous color map

#Store the concentrations for each integer position in a matrix
m = 4000
conc_matrix = np.zeros((m, m)) #we will display from [-1000, -1000] to [3000, 3000]
for i in range(m):
    for j in range(m):
        conc_matrix[i][j] = math.log(calculate_conc([i - 1000, j - 1000]))


#Simulate the gradient distribution, plot as a heatmap
ax.imshow(conc_matrix.T, cmap=cmap_color, interpolation='nearest', extent = [plot_lim_x[0], plot_lim_x[1], plot_lim_y[0], plot_lim_y[1]], origin = 'lower')

ax.set_title("Biased walk with chemoattractant".format(), x = 0.5, y = 0.87)
ax.set_xlim(plot_lim_x[0], plot_lim_x[1])
ax.set_ylim(plot_lim_y[0], plot_lim_y[1])
ax.set_xlabel("poisiton in μm")
ax.set_ylabel("poisiton in μm")

#ABP_step_biased_linear()
fig.tight_layout()
plt.show()