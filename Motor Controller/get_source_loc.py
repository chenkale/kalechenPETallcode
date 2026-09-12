import cv2
import numpy as np
import matplotlib.pyplot as plt
import time
import os

script_dir = os.path.dirname(__file__)

def fit_image_to_canvas(image_path, canvas_size=255, save_path= script_dir + '/resized_img.png'):
    # Load image
    img = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
    if img is None:
        raise ValueError("Image not found or unreadable.")

    # Convert image to gray scale if not
    if len(img.shape) == 3:
        img_gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    else:
        img_gray = img

    # Threshold to get binary shape (invert so black shapes are detected)
    _, thresh = cv2.threshold(img_gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)

    # Find shapes in image
    contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        raise ValueError("No shape detected.")

    # Use largest shape
    largest_contour = max(contours, key=cv2.contourArea)
    x, y, w, h = cv2.boundingRect(largest_contour)

    # Crop image around shape
    cropped_shape = img_gray[y:y+h, x:x+w]

    # Resize shape to fit canvas while maintaining ratios
    scale = min(canvas_size / w, canvas_size / h)
    new_w = int(w * scale)
    new_h = int(h * scale)
    resized_shape = cv2.resize(cropped_shape, (new_w, new_h), interpolation=cv2.INTER_AREA)

    # Create white canvas
    canvas = np.ones((canvas_size, canvas_size), dtype=np.uint8) * 255

    # Center resized shape on canvas
    x_offset = (canvas_size - new_w) // 2
    y_offset = (canvas_size - new_h) // 2
    canvas[y_offset:y_offset + new_h, x_offset:x_offset + new_w] = resized_shape

    # Save and return
    cv2.imwrite(save_path, canvas)
    return canvas

# Orders points by proximity to one another
def sortbydis(poslist, option = 0):
    sortedlist = [[0,0]]
    for i in range(len(poslist)):
        shortdis = 1000000000000
        temppos = 0
        for j in range(len(poslist)): # Checks distance between every point and the last point in the current list
            distance = np.sqrt((sortedlist[-1][1]-poslist[j][1])**2 + (sortedlist[-1][0]-poslist[j][0])**2) # Distance
            if(distance < shortdis): 
                shortdis = distance
                temppos = j
        sortedlist.append(poslist[temppos]) # Creates new sorted list
        poslist.remove(poslist[temppos])
    
    del sortedlist[0]
    
    if (option == 1):
        for i in range(len(sortedlist)):
            sortedlist[i][0] = sortedlist[i][0] - 125 
            sortedlist[i][1] = sortedlist[i][1] - 125
    
    return sortedlist

# Restricts number of points to user specified amount
def limit_points(num, old_list):

    while (len(old_list) < num):
        old_list = addPoints(old_list)

    new_list = []
    num_remove = len(old_list) - num
    remove_array = np.linspace(0, len(old_list) - 1, num_remove, dtype = int)


    for i, element in enumerate(old_list):
        if not (i in remove_array):
            new_list.append(element)
            
#    print("Len new list: " + str(len(new_list)))
    
    return new_list
            
            
        
    
def num_points(data_points, num):
    len(data_points)

def store_list(data_list):
        with open(script_dir + "/source_locations", 'w') as file:
            for row in data_list:
                row_str = ' '.join(map(str, row))
                file.write(row_str + '\n')
                

def totDistance(points):
    sumDis = 0
    for i in range(1, len(points)):
        distance = np.sqrt((points[i][1]-points[i-1][1])**2 + (points[i][0]-points[i-1][0])**2)
        sumDis = sumDis + distance
    return sumDis

def addPoints(points):
    midPoints = [[0,0]]
    for i in range(len(points)-1,0,-1):
        x_coord = (points[i][0] + points[i-1][0])/2
        y_coord = (points[i][1] + points[i-1][1])/2
        midPoints.append([x_coord,y_coord])
    
    del midPoints[0]

    newPoints = [[0,0]]
    
    for i, point in enumerate(points):
        #print("i is: " + str(i))
        newPoints.append(point)
        if (i != len(points) - 1):
            newPoints.append(midPoints[i])
        
        
    del newPoints[0]
    
    newPoints = sortbydis(newPoints)
    
#    print("Len new points: " + str(len(newPoints)))
    
    return newPoints
    

def main():
        
        

        
    #Gets picture file
    files = os.listdir(script_dir + '/image')
    if (len(files) > 0):
        if ('png' in files[0]):
            file_name = files[0]
        else:
            print("Please upload png")
            return
    else:
        print("No image file found/usable")
        return

    #Gets dimensions of shape
    user_scale_num = 300
    while (user_scale_num > 1 or user_scale_num <= 0):
        user_scale_num_str = input("Enter scale of shape (0-1): ")
        try:
            user_scale_num = float(user_scale_num_str)
            if (user_scale_num > 1 or user_scale_num <= 0):
                print("Please enter a number between 0 and 1")
        except ValueError:
            print("Please enter an number")
           


    fit_image_to_canvas(script_dir + '/image/' + file_name)
    img = cv2.imread(script_dir + "/resized_img.png")
    img = cv2.flip(img, 0)
    gray_img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    height, width = gray_img.shape
    poslist = []

    _, binary = cv2.threshold(gray_img, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
    contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)

    corners = contours[0].squeeze()
    

    
    for c in corners:
        x, y = c.ravel()
        poslist.append([x,y])
        img = cv2.circle(img, center=(x, y), radius=1, 
                        color=(0, 0, 255), thickness=-1)
        
        
    
    print("Suggested Points: " + str(len(poslist)))

    sorted_points = sortbydis(poslist)
    
    
                        
    #Gets amount of coordinates motor will use from user (depends on prefered resolution)
    #user_num_points = len(sorted_points) + 1
    user_num_points = 10000000000000
    while (len(sorted_points) - user_num_points < -10000000):
        user_num_points_str = input("Enter number of points (Press enter to not specify): ")
        try:
            user_num_points = int(user_num_points_str)
        except ValueError:
            user_num_points = 0

    
    if (user_num_points != 0):
        sorted_points = limit_points(user_num_points, sorted_points)
        
    sorted_points = sortbydis(sorted_points, 1)
    sorted_points = [[elem * user_scale_num for elem in row] for row in sorted_points]
        
    # Writes points to file 
    store_list(sorted_points)


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
        for i in range(len(sorted_points)): 
                x_data.append(sorted_points[i][0])
                y_data.append(sorted_points[i][1]) 

                ax.clear() # Clear previous plot
                ax.plot(x_data, y_data, 'ro') # Plot all points as red circles
                ax.set_xlim(-width/2, width/2) # Set x-axis limits
                ax.set_ylim(-height/2, height/2) # Set y-axis limits

                fig.canvas.draw()
                fig.canvas.flush_events()
                time.sleep(0.00005) # Pause time (seconds)
                
                
                
#    for p in sorted_points:
#        x, y = p
#        img = cv2.circle(img, center=(x, y), radius=1, 
#                        color=(0, 0, 255), thickness=-1)
                        
    cv2.imwrite(script_dir + "/traced_img.png", img)

    # Gets total distance between points
    distance = totDistance(sorted_points)

    #print("distance: " + str(distance))
    print("Num Points: " + str(len(sorted_points)))


    return distance
    
#main()

