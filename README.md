# Simulation for slime mold growing in response to food sources
This GitHub repository contains the scripts necessary to simulate slime mold moving in response to food by using an Active Brownian Particle model to simulate their search and introducing a chemoattractant force.

## Contributors
Jesse Rodriguez Reyes, Justin Couvrette, Salmanali Mohammad, Tyler Couvrette

## Getting Started
The following packages were used when running the code:
python 3.8.8
numpy version 1.20.1
matplotlib 3.3.4
scipy 1.6.2
math 

### Before you begin:
Make sure all mentioned packages are installed. Run the first code chunk to load all packages mentioned. There are also two markdown files included. One is BIOE234_final, which is a fully functional markdown that includes all objectives we successfully have running so far. The second is BIOE234_final_work_in_progress, which has objectives that were/are being acively worked on and includes detailed explanations on what the objective is and how we are trying to ahcieve it. The work in progress objectives currently do not work. 

## Utilizing the script
The simulations are run using a single markdown file. The markdown file is segmented into multiple different sections, with each section achieving a different goal. The sections, called objectives, are listed below.

### Objective 1: Model random walkers trajectories that move in response to a single chemoattractant
No input needed initially.
Running file will produce 10 completely random walkers in one section and 10 random walkers that are drawn to an attractive force in the next.

The paramerters of the ABP model can be modified to change the behavior of the random walk aspect. These include the velocity of the walkers, the timesteps used, and the rotational and translational diffusion. The number of walkers and steps can be changed at the beginning as well. The nutrients location can be altered by changing the nutrient_position. The strength of the chemoattractant can be altered by changing the nutrient exponent of start exponent. Increasing the nutrient exponent will increase the concentration of the chemoattractant at the nutrient source, decreasing it will do the opposite. Increasing the start exponent will increase the concentration of the chemoattractant at the start of the walk, decreasing it will do the opposite. These can be changed to modify the strength of the chemoattractant force through the walk. Changing x_scaling will control where the random walkers are initialized, from 0 to 1*x_scaling. 

### Objective 2: Model random walkers trajectories that move in response to chemoattractants scattered along a straight line
No input needed
Running file will produce 10 random walkers that are drawn to an uniformly upwards attractive force where there is no x component to the chemotrractant.

The paramerters of the ABP model can be modified to change the behavior of the random walk aspect. These include the velocity of the walkers, the timesteps used, and the rotational and translational diffusion. The number of walkers and steps can be changed at the beginning as well. The nutrients location can be altered by changing the nutrient_position. The strength of the chemoattractant can be altered by changing the nutrient exponent of start exponent. Increasing the nutrient exponent will increase the concentration of the chemoattractant at the nutrient source, decreasing it will do the opposite. Increasing the start exponent will increase the concentration of the chemoattractant at the start of the walk, decreasing it will do the opposite. These can be changed to modify the strength of the chemoattractant force through the walk. Changing x_scaling will control where the random walkers are initialized, from 0 to 1*x_scaling. 


