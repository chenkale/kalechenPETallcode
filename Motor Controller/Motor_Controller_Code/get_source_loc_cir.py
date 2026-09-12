import numpy as np
import matplotlib.pyplot as plt
import time
import os

starting_angle = np.pi / 2


script_dir = os.path.dirname(__file__)

# Orders points by proximity to one another


# Restricts number of points to user specified amount
  
    

def store_list(data_list):
        with open(script_dir + "/source_locations", 'w') as file:
            for row in data_list:
                row_str = ' '.join(map(str, row))
                file.write(row_str + '\n')
                

def totDistance(points):
    sumDis = 0
    for i in range(1, len(points)):
        distance = np.sqrt((points[i][1]-points[i-1][1])**2 + (points[i][0]-points[i-1][0])**2)
        if (distance < 10):
            sumDis = sumDis + distance
    return sumDis


def main(radius):
        
        #Gets dimensions of shape
    user_radius_num = -1
    '''
    while (user_scale_num < 0):
        user_radius_num_str = input("Enter a radius (mm): ")
        try:
            user_radius_num = float(user_radius_num_str)
        except ValueError:
            print("Please enter an number")    
    '''
    
    thetas = np.linspace(starting_angle,2*np.pi + starting_angle,2000)
    i = 0
    #for theta in thetas:
    #    if (-1.12011132679 > theta > -2.02148132679 or -4.26170398 > theta > -5.16307398):
    #        thetas = np.delete(thetas, i)
    #         i = i-1
    #    i = i + 1
    circle_list =[[0,0]]
    for theta in thetas:
        circle_list.append([radius*np.cos(theta),radius*np.sin(theta)])    
    del circle_list[0]
#    print(circle_list)
    
    # Writes points to file 
    store_list(circle_list)


    # Checks whether user wants to plot points
    plot_bool_str = input("Plot points yes (1) or no (0): ")
    if (plot_bool_str == '1'):
        plot_bool = True
    else:
        plot_bool = False

        
    # Plots all points if requested (intented to double check becfore starting a run)

    if (plot_bool):

        plt.ion()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        x_data = []
        y_data = []
        for i in range(len(circle_list)): 
                x_data.append(circle_list[i][0])
                y_data.append(circle_list[i][1]) 

                ax.clear() # Clear previous plot
                ax.plot(x_data, y_data, 'ro') # Plot all points as red circles
                ax.set_xlim(-400/2, 400/2) # Set x-axis limits
                ax.set_ylim(-400/2, 400/2) # Set y-axis limits

                fig.canvas.draw()
                fig.canvas.flush_events()
                time.sleep(0.00005) # Pause time (seconds)
                
                
                
#    for p in circle_list:
#        x, y = p
#        img = cv2.circle(img, center=(x, y), radius=1, 
#                        color=(0, 0, 255), thickness=-1)
                        

    # Gets total distance between points
    distance = totDistance(circle_list)

    #print("distance: " + str(distance))
    print("Num Points: " + str(len(circle_list)))


    return distance

