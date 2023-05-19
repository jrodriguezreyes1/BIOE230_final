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
