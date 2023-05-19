import numpy as np
import matplotlib.pyplot as plt
from numpy.random import random as rand
from scipy import stats


import math
from matplotlib import colors# Objective 1: Model random walkers trajectories that move in response to a single chemoattractant# Initialize parameters of ABP model and locations first

#ABP model parameters: ballistic velocity, time step, rotational diffusion constant, translational diffusion constant
vel =5.0; dt = 1; Drot = 0.1; Dtrans = 0.1;
# vel =5.0; dt = 1; Drot = 0.1; Dtrans = 0.1; #These work pretty well, will see if other work better

# initialize arrays that store x,y and theta values, 
#as well as initial particle position and angle
num_steps = 500;
num_walks = 10;
x_scaling = 3000;
 
nutrient_center = [1500, 1500] # point where nutrient is located
nutrient_exponent, start_exponent = 6, 2 # Initial concentration and 10^3 uM
proximity = 0 # distance from start to center, will be calculated later# Define completely random walkers where steps are taken using ABP

def ABP_step_rand():
    for j  in range(num_walks):
        
        xvec=np.zeros(0); yvec=np.zeros(0); thetavec = np.zeros(0);
        x=x_scaling*(rand()); y = 0.0; theta = (2*np.pi)*rand();
        
        # inner for loop is each step for a given walker/trajectory
        for i in range(num_steps):
                # calculate diffusive/random steps. For the x- and y-,we generate 
                #a random number between -1 & 1 to ensure that the walker can step in both directions(up/down and left/right).

                dx = np.sqrt(2*Dtrans*dt)*2*(rand(1)-0.5); 
                dy= np.sqrt(2*Dtrans*dt)*2*(rand(1)-0.5); 
                dtheta = np.sqrt(2*Drot*dt)*2*(rand(1)-0.5);

                # update coordinates (including ballistic step)
                x += vel*dt*np.cos(theta) + dx;
                y += vel*dt*np.sin(theta) + dy;


                # store successive positions in arrays
                xvec = np.append(xvec,x); yvec = np.append(yvec,y); 
                # update the angle and store in array
                theta += dtheta;
                thetavec = np.append(thetavec, theta);
        plt.plot(xvec,yvec, '-');
        plt.plot(xvec[-1],yvec[-1],'ro', markersize = 8)
        plt.plot(xvec[0],yvec[0],'bo', markersize = 8)
    
#ABP_step_rand()# Define biased random walker with chemoattractant where steps are taken using ABP


# Calculates the distance from the random walkers current position (a) to the focal point of the nutrient (b)
def distance_to_nutrient(a,b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)


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

    x_force = (xnutr - x)*concentration*10**-6; # as concentration increases, make walker move towards 
    y_force = (ynutr - y)*concentration*10**-6; # as concentration increases, make walker move towards 
    
    x_combined = x+x_force;
    y_combined = y+y_force;
    return x_combined, y_combined, x_force, y_force, proximity, concentration


def ABP_step_biased():

    for j  in range(num_walks):
        xvec=np.zeros(0); yvec=np.zeros(0); thetavec = np.zeros(0);
        x=x_scaling*(rand()); y = 0.0; theta = (2*np.pi)*rand();
        
        proximity = distance_to_nutrient([x,y],nutrient_center); # initialize proximity to nutrient at start
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
                #dtheta = np.sqrt(2*Drot*dt)*2*(rand(1)-0.5);

                # update coordinates (including ballistic step)
                x += vel*dt*np.cos(theta) + dx;
                y += vel*dt*np.sin(theta) + dy;

                # store successive positions in arrays
                xvec = np.append(xvec,x); yvec = np.append(yvec,y); 
                # update the angle and store in array
                theta += dtheta;
                thetavec = np.append(thetavec, theta);
                
                
                pos_comb = gradient_force(x,y,nutrient_center[0],nutrient_center[1])
    
                # store successive positions in arrays
                x = pos_comb[0];
                y = pos_comb[1];
                                
        plt.plot(xvec,yvec, '-');
        plt.plot(xvec[-1],yvec[-1],'ro', markersize = 8)
        plt.plot(xvec[0],yvec[0],'bo', markersize = 8)
        plt.plot(nutrient_center[0],nutrient_center[1],'g*', markersize = 16)

#ABP_step_biased()# Make plot of chemoattractant gradient with random walk
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

ax.plot(nutrient_center[0], nutrient_center[1], 'bX', markersize = 8) #Mark the highest concentration point [1500, 1500]
ax.set_title(" Completely Random walk".format(), x = 0.5, y = 0.87)
ax.set_xlim(plot_lim_x[0], plot_lim_x[1])
ax.set_ylim(plot_lim_y[0], plot_lim_y[1])
ax.set_xlabel("poisiton in μm")
ax.set_ylabel("poisiton in μm")

#ABP_step_rand()
fig.tight_layout()
plt.show()# Make plot of chemoattractant gradient with biased random walk that responds to chemoattractant
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

ax.plot(nutrient_center[0], nutrient_center[1], 'bX', markersize = 8) #Mark the highest concentration point [1500, 1500]
ax.set_title("Biased walk with chemoattractant".format(), x = 0.5, y = 0.87)
ax.set_xlim(plot_lim_x[0], plot_lim_x[1])
ax.set_ylim(plot_lim_y[0], plot_lim_y[1])
ax.set_xlabel("poisiton in μm")
ax.set_ylabel("poisiton in μm")

#ABP_step_biased()
fig.tight_layout()
plt.show()#  Objective 2.1: Model random walkers trajectories that move in response to chemoattractants scattered along a straight line# ABP model parameters: ballistic velocity, time step, rotational diffusion constant, translational diffusion constant

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
    return x_combined, y_combined, x_force, y_force, proximity, concentration# Define biased random walker with chemoattractant where steps are taken using ABP and chemoattractant is set across
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

#ABP_step_biased_linear()# Make plot of chemoattractant gradient with biased random walk that responds to chemoattractant
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
plt.show()# Make plot of chemoattractant gradient with random walk
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

ax.set_title(" Completely Random walk".format(), x = 0.5, y = 0.87)
ax.set_xlim(plot_lim_x[0], plot_lim_x[1])
ax.set_ylim(plot_lim_y[0], plot_lim_y[1])
ax.set_xlabel("poisiton in μm")
ax.set_ylabel("poisiton in μm")

#ABP_step_rand()
fig.tight_layout()
plt.show()# Objective 2.2: Model random walkers trajectories that move in response to chemoattractants scattered along a line defined by an equation

### The objective here is to be able to define a chemoattractant line using some function and have the walkers move in response to the proximity to the chemoattractant, with both a x and y attractive force instead of just the y attratcive force in 2.1. We need to add a function that calculates the distance to the nearest point on the chemoattractant line then used the exisiting gradient_force function used in objective 1 to calculate the x and y chemoattractant force towards that nearest point. The purpose of this is for materials design applications, where the walkers will move towards the lower region of a parabola for example which can be used to design support networks for an actual structure/material that is the same as that parabola. An example of this is seen near the end of the powerpoint in this GitHub

### The code is currently nonfunctional and a work in progress. More work will be done on it beyond the BIOE230 course.# ABP model parameters: ballistic velocity, time step, rotational diffusion constant, translational diffusion constant

# Now considers distance to nearest point on line representing th chemoattractant, which isn't necessarily straight 




# Initialize parameters of ABP model and locations first

# ABP model parameters: ballistic velocity, time step, rotational diffusion constant, translational diffusion constant
vel =5.0; dt = 1; Drot = 0.1; Dtrans = 0.1;

# initialize arrays that store x,y and theta values, 
# as well as initial particle position and angle
num_steps = 1000;
num_walks = 10;
x_scaling = 2000;

nutrient_center = [0, 1500] # point where nutrient is located
x_nutrient = np.arange(-10, 10, 2, dtype=int)
print(x_nutrient)
y_nutrient = 1.5*(x_nutrient);
plt.plot(x_nutrient,y_nutrient)
nutrient_center = [x_nutrient,y_nutrient]

nutrient_exponent, start_exponent = 6, 2.5 # Concentration is 10^exponent for each respective one
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
    return x_combined, y_combined, x_force, y_force, proximity, concentration# Define biased random walker with chemoattractant where steps are taken using ABP and chemoattractant is set across
# a line at the top



# ABP model parameters: ballistic velocity, time step, rotational diffusion constant, translational diffusion constant
vel =5.0; dt = 1; Drot = 0.1; Dtrans = 0.1;

# initialize arrays that store x,y and theta values, 
# as well as initial particle position and angle
num_steps = 1000;
num_walks = 10;
x_scaling = 2000;


nutrient_center = [0, 1500] # point where nutrient is located
x_nutrient = np.arange(-10, 10, 2, dtype=int)
print(x_nutrient)
y_nutrient = 1.5*(x_nutrient);
plt.plot(x_nutrient,y_nutrient)
nutrient_center = [x_nutrient,y_nutrient]
print(nutrient_center)


nutrient_exponent, start_exponent = 6, 2.5 # Concentration is 10^exponent for each respective one
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
                #dtheta = np.sqrt(2*Drot*dt)*2*(rand(1)-0.5);

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

#ABP_step_biased_linear()# Objective 3: Build random walkers simultaneously to later enable attraction and repulsion towards other walkers as they move

## The goal here is to invert the loop order by making the walkers each take a step at a time before moving on to the next walker to do one step so they get walk together, just like slime mold will explore with multiple walkers moving at the same time as opposed to one full walk being taken before the next full walk is made. This is to later add attractive and repulsive forces to the walkers as they are made based on proximity to the existing walker trajectories. This force will be repulsive early on and attractive later on to more closely mimic the behavior of slime mold. 

### The code is currently a work in progress and nonfunctional. More work will be done on it beyond the BIOE230 course. # ABP model parameters: ballistic velocity, time step, rotational diffusion constant, translational diffusion constant

# Removed x term from distance calculation so only y vector of chemoattractant is considered



# Initialize parameters of ABP model and locations first

# ABP model parameters: ballistic velocity, time step, rotational diffusion constant, translational diffusion constant
vel =8.0; dt = 1; Drot = 0.03; Dtrans = 3;
#vel =5.0; dt = 1; Drot = 0.1; Dtrans = 0.1; #original parameters that work well in case the new ones get screwed up

# initialize arrays that store x,y and theta values, 
# as well as initial particle position and angle
num_steps = 100;
num_walks = 5;
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
    return x_combined, y_combined, x_force, y_force, proximity, concentration# Define biased random walker with chemoattractant where steps are taken using ABP and chemoattractant is set across
# a line at the top


def ABP_step_biased_linear():

    xvec=np.zeros(0); yvec=np.zeros(0); thetavec = np.zeros(0);
    x=x_scaling*(rand()); y = 0.0; theta = (2*np.pi)*rand();

    proximity = distance_to_nutrient([x,y],nutrient_center); #initialize proximity to nutrient at start
#         print('Initial proximity =',proximity)
    concentration = calculate_conc([x,y]) # initialize concentration of nutrient at start
#         print('Initial concentration =',concentration,'\n')

    # inner for loop is each step for a given walker/trajectory        
    for i in range(num_steps):
            
        for j  in range(num_walks):
            
            print("step = ", i, "\n walk =", j)
            if i == 0:
                xvec=np.zeros(0); yvec=np.zeros(0); thetavec = np.zeros(0);
                x=x_scaling*(rand()); y = 0.0; theta = (2*np.pi)*rand();
                
                proximity = distance_to_nutrient([x,y],nutrient_center); #initialize proximity to nutrient at start
            #         print('Initial proximity =',proximity)
                concentration = calculate_conc([x,y]) # initialize concentration of nutrient at start
            #         print('Initial concentration =',concentration,'\n')

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
                
            else:

                proximity = distance_to_nutrient([x,y],nutrient_center); #initialize proximity to nutrient at start
            #         print('Initial proximity =',proximity)
                concentration = calculate_conc([x,y]) # initialize concentration of nutrient at start
            #         print('Initial concentration =',concentration,'\n')

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
                
                #return x, y
                

    plt.plot(xvec,yvec, '-');
    plt.plot(xvec[-1],yvec[-1],'ro', markersize = 8)
    plt.plot(xvec[0],yvec[0],'bo', markersize = 8)
#         plt.plot(nutrient_center[0],nutrient_center[1],'g*', markersize = 16)

#ABP_step_biased_linear()### BACKUP CODE IN CASE THE CODE ABOVE GETS IRREPARABLY BROKEN# Define biased random walker with chemoattractant where steps are taken using ABP and chemoattractant is set across
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

#ABP_step_biased_linear()