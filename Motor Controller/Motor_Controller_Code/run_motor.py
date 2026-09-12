# In order for this code to run properly please manually re-position the arm of the motor
# so that it is aligned with the corner closest to the power source.
# Also note to run this code please install the pyserial module : pip install pyserial

import tkinter as tk
from tkinter.messagebox import askyesno
import platform
import serial
import serial.tools.list_ports
from serial.tools import list_ports
import numpy as np
import os
import time
import get_source_loc_cir
from datetime import datetime

centerx, centery = 197.5, 215

def wait_for_text(expected_text):
    entered_text = ''
   
    while entered_text != expected_text:
        # Use Javascript to prompt for text input in Jupyter Notebook
        print('Stopped Motor')
        print('Type \'go\' and press enter to restart motor')
        entered_text = input("Enter text: ")
    print('Restarting Motor')


# enter the time you need the source to stop in any position and the array of coordinates used for the run
def getTime(timeNeeded, array):
    return timeNeeded/ len(array)

def getFeedRate(timeNeeded,distance):
    return (distance/timeNeeded)*60

def getDistance(i, points):
    return np.sqrt((points[i][1]-points[i-1][1])**2 + (points[i][0]-points[i-1][0])**2)

# Searches through usb devices on computer to locate Aufero Laser Engraver 2
def find_aufero_port(baud=115200, timeout=.3):
    
    # Specifies the type of port's to check
    TARGET_VID_PID = {
    (0x1A86, 0x7523),  # CH340
    (0x10C4, 0xEA60),  # CP2102
    }

    for port in list_ports.comports():
        # Only check ports with known VID/PID
        desc = port.description or ""
        device = port.device or ""

        # Skip Bluetooth ports
        if "Bluetooth" in desc or "bluetooth" in desc or "rfcomm" in device or "cu.Bluetooth" in device:
            continue
            
        # Skips port types not specified    
        #if (port.vid, port.pid) not in TARGET_VID_PID:
        #    continue

        try:
            with serial.Serial(port.device, baud, timeout=timeout) as ser: 
                ser.reset_input_buffer()
                ser.write(b"\r\n")  # Safe wake-up command
                ser.flush()
                ser.write(b"?\n")  # GRBL status check
                response = ser.read(100).decode(errors='ignore')

                if any(keyword in response for keyword in ["Grbl", "<Idle", "<Run", "ok"]): # Possible responses from Aufero
                    print(f"Aufero detected on {port.device}")
                    return port.device
                else:
                    print(f"No GRBL response from {port.device}")

        except (serial.SerialException, OSError) as e:
            print(f"Could not open {port.device}: {e}")
            continue

    print("Aufero Laser 2 not found.")
    return None

port = find_aufero_port()
#print("Port Found: " + str(port))
#print("Port Name: " + str(port.device))

# logistical code for connecting to the driver
platform.system() =='Linux'
ser = serial.Serial()
ser.port = str(port)
ser.baudrate = 115200              # speed of communication between computer and port
ser.bytesize = serial.EIGHTBITS    # number of bits per bytes
ser.parity = serial.PARITY_NONE    # parity checks wheter or not data is lost or overwritten... we dont need
ser.stopbits = serial.STOPBITS_ONE # number of stop bits
ser.timeout = 1                    # time til time out
ser.xonxoff = False                # disable software flow control
ser.rtscts = False                 # disable hardware (RTS/CTS) flow control
ser.dsrdtr = False                 # disable hardware (DSR/DTR) flow control
ser.writeTimeout = 2               # timeout for write
script_dir = os.path.dirname(__file__)

def shapeMotion(points, timeOfRun):
    ser.close()
    ser.open()
    ser.write(b"\r\n\r\n")
    time.sleep(2)     # Wait for Printrbot to initialize
    ser.flushInput()  # Flush startup text in serial inp

    string = "\r\n\r\n"
    ser.write(str.encode(string))
    time.sleep(2)     # Wait for Printrbot to initialize
    ser.flushInput()  # Flush startup text in serial inp
    ### start up 
    string = 'M18 ; \n'
    ser.write(str.encode(string))
    time.sleep(2)     # in units of seconds 

    ### set absolute positioning ###
    string = 'G90 ; \n'
    ser.write(str.encode(string))
    time.sleep(2)

    ### set units to mm ###
    string = 'G21 ; \n' 
    ser.write(str.encode(string))
    time.sleep(2)

    times = getTime(timeOfRun,points)
    i = 0
    
    User_Input = 1
    while (User_Input == 1):

        for point in points:
            if i ==0: # Positions motor for run
                time.sleep(2)
                string = 'G17 X{} Y{} F{} ; \n'.format(198,215,5000)
                ser.write(str.encode(string))
                string = 'G17 X{} Y{} F{} ; \n'.format(198+point[0],215+point[1],5000)
                ser.write(str.encode(string))
                wait_for_text('go')
                i += 1
            else:
            # Loops through points on shape
            # Note that this center point is setup specific and was found by trial-and-error
            # This center point places the center of the circle at the center of the scanner when
            # the 2D stage is mounted onto the top of the scanner 
                string = 'G17 X{} Y{} F100 ; \n'.format(198+point[0],215+point[1])
                ser.write(str.encode(string))
                time.sleep(times)
        
        # Checks whether user wants to run again
        user_input_str = '2'
        while (user_input_str != '1' and user_input_str != '0'):
            user_input_str = input("Run Again (1) or Stop (0): ")
            try:
                User_Input = int(user_input_str)
                if (User_Input != 1 and User_Input != 0):
                    print("Enter 1 or 0")
            except ValueError:
                print("Enter an integer")
        i = 1
        
    # shift back to origin (0,0) position after text prompt
    wait_for_text('go')
    string = 'G17 X{} Y{} F5000 ; \n'.format(0,0)
    ser.write(str.encode(string))
    time.sleep(2)

def linearShapeMotion(points, timeOfRun, distance, lineRun = False):
    ser.close()
    ser.open()
    ser.write(b"\r\n\r\n")
    time.sleep(2)     # Wait for Printrbot to initialize
    ser.flushInput()  # Flush startup text in serial inp

    string = "\r\n\r\n"
    ser.write(str.encode(string))
    time.sleep(2)     # Wait for Printrbot to initialize
    ser.flushInput()  # Flush startup text in serial inp
    ### start up 
    string = 'M18 ; \n'
    ser.write(str.encode(string))
    time.sleep(2)     # in units of seconds 

    ### set absolute positioning ###
    string = 'G90 ; \n'
    ser.write(str.encode(string))
    time.sleep(2)

    ### set units to mm ###
    string = 'G21 ; \n' 
    ser.write(str.encode(string))
    time.sleep(2)
    
    feedrate = getFeedRate(timeOfRun,distance)
    print("Feedrate: " + str(feedrate))
    print("Distance: " + str(distance) + " mm")
    print("Time Calculated (From Feedrate): " + str(distance/(feedrate/60)))

    i = 0

    User_Input = 1
    while (User_Input == 1):
        for point in points:

            if i ==0:
                time.sleep(2)
                string = 'G17 X{} Y{} F{} ; \n'.format(centerx,centery,5000)
                ser.write(str.encode(string))
                string = 'G17 X{} Y{} F{} ; \n'.format(centerx+point[0],centery+point[1],5000)
                ser.write(str.encode(string))
                ser.write(b"M400\n")
                time.sleep(0.1)

                wait_for_text('go')
                i += 1
            else:
                string = 'G1 X{} Y{} F{} ; \n'.format(centerx+point[0],centery+point[1],feedrate)
                
                ser.write(str.encode(string))

                if(i%20 == 0):
                    time.sleep(((timeOfRun)/len(points))*18)
                i+=1
        if(lineRun):
            for i in range(0, int(timeOfRun/100)):
                time.sleep(99)
                ser.write(b"G4 P0\n")
                
        string = 'G1 X{} Y{} F{} ; \n'.format(centerx+points[0][0],centery+points[0][1],5000)
        ser.write(str.encode(string))
        
        user_input_str = '2'
        while (user_input_str != '1' and user_input_str != '0'):
            user_input_str = input("Run Again (1) or Stop (0): ")
            try:
                User_Input = int(user_input_str)
                if (User_Input != 1 and User_Input != 0):
                    print("Enter 1 or 0")
            except ValueError:
                print("Enter an integer")
        i = 1
                
        

    wait_for_text('go')
    string = 'G17 X{} Y{} F5000 ; \n'.format(0,0)
    ser.write(str.encode(string))
    time.sleep(2)

def discreteLinearShapeMotion(start_point, end_point, num_points, timeOfRun, lineRun=False):
    """
    Moves quickly between discrete points in stop-and-go fashion. 
    Spends most of the time stopped at each point, dividing time equally.
    """
    ser.close()
    ser.open()
    ser.write(b"\r\n\r\n")
    time.sleep(2)
    ser.flushInput()

    string = "\r\n\r\n"
    ser.write(str.encode(string))
    time.sleep(2)
    ser.flushInput()

    # Start up and config
    ser.write(str.encode('M18 ; \n'))
    time.sleep(2)
    ser.write(str.encode('G90 ; \n'))
    time.sleep(2)
    ser.write(str.encode('G21 ; \n'))
    time.sleep(2)

    if num_points < 2:
        raise ValueError("num_points must be at least 2")
    x0, y0 = start_point
    x1, y1 = end_point
    xs = np.linspace(x0, x1, num=num_points)
    ys = np.linspace(y0, y1, num=num_points)
    points = [(x, y) for x, y in zip(xs, ys)]

    # Calculate fast feedrate for short jump between points
    move_feedrate = 5000  # very high for fast movement, adjust as safe
    dwell_time = float(timeOfRun) / num_points  # time spent at each point (sec)

    print("Stop-and-go: {} discrete points, {:.2f}s/point, fast feedrate: {}".format(num_points, dwell_time, move_feedrate))

    User_Input = 1
    while User_Input == 1:
        # Move to origin (center)
        ser.write(str.encode('G17 X{} Y{} F{} ; \n'.format(centerx, centery, move_feedrate)))
        time.sleep(0.2)
        # Optionally wait for user to start
        wait_for_text('go')

        for idx, point in enumerate(points):
            # Move quickly to next point
            ser.write(str.encode('G1 X{} Y{} F{} ; \n'.format(centerx + point[0], centery + point[1], move_feedrate)))
            ser.write(b"M400\n")
            time.sleep(0.1)  # Give motion command a moment to start
            # Dwell (stop) at location
            time.sleep(dwell_time)
        
        if lineRun:
            for _ in range(0, int(timeOfRun / 100)):
                time.sleep(99)
                ser.write(b"G4 P0\n")

        # Return quickly to 1st point
        ser.write(str.encode('G1 X{} Y{} F{} ; \n'.format(centerx + points[0][0], centery + points[0][1], move_feedrate)))
        time.sleep(0.2)

        user_input_str = '2'
        while (user_input_str != '1' and user_input_str != '0'):
            user_input_str = input("Run Again (1) or Stop (0): ")
            try:
                User_Input = int(user_input_str)
                if (User_Input != 1 and User_Input != 0):
                    print("Enter 1 or 0")
            except ValueError:
                print("Enter an integer")

    wait_for_text('go')
    string = 'G17 X{} Y{} F5000 ; \n'.format(0, 0)
    ser.write(str.encode(string))
    time.sleep(2)

def discretePoints(points, user_time_num):
    """
    Moves the motor through a list of (x, y) points, pausing at each for dwell_time seconds.
    Runs similar to discreteLinearShapeMotion, but takes all points from user.
    """

    ser.close()
    ser.open()
    ser.write(b"\r\n\r\n")
    time.sleep(2)
    ser.flushInput()

    string = "\r\n\r\n"
    ser.write(str.encode(string))
    time.sleep(2)
    ser.flushInput()

    # Start up and config
    ser.write(str.encode('M18 ; \n'))
    time.sleep(2)
    ser.write(str.encode('G90 ; \n'))
    time.sleep(2)
    ser.write(str.encode('G21 ; \n'))
    time.sleep(2)

    move_feedrate = 5000  # very high for fast movement, adjust as safe
    dwell_time = user_time_num / len(points)
    print("Discrete points: {} total, {:.2f}s per point, feedrate: {}".format(
        len(points), dwell_time, move_feedrate))

    User_Input = 1
    while User_Input == 1:
        # Move to origin (center)
        ser.write(str.encode('G17 X{} Y{} F{} ; \n'.format(centerx, centery, move_feedrate)))
        time.sleep(0.2)
        # Optionally wait for user to start
        wait_for_text('go')

        for idx, point in enumerate(points):
            # Move quickly to next point
            ser.write(str.encode('G1 X{} Y{} F{} ; \n'.format(centerx + point[0], centery + point[1], move_feedrate)))
            ser.write(b"M400\n")
            time.sleep(0.1)  # Give motion command a moment to start
            # Dwell (stop) at location
            time.sleep(dwell_time)
        
        # Return quickly to 1st point
        if points:
            ser.write(str.encode('G1 X{} Y{} F{} ; \n'.format(centerx + points[0][0], centery + points[0][1], move_feedrate)))
            time.sleep(0.2)

        user_input_str = '2'
        while (user_input_str != '1' and user_input_str != '0'):
            user_input_str = input("Run Again (1) or Stop (0): ")
            try:
                User_Input = int(user_input_str)
                if (User_Input != 1 and User_Input != 0):
                    print("Enter 1 or 0")
            except ValueError:
                print("Enter an integer")

    wait_for_text('go')
    string = 'G17 X{} Y{} F5000 ; \n'.format(0, 0)
    ser.write(str.encode(string))
    time.sleep(2)


def read_loc():
    data_list = []
    with open(script_dir + "/source_locations", 'r') as file:
        for line in file:
            row = list(map(float, line.strip().split()))
            data_list.append(row)
    return data_list
    
def circularMotion(points, timeOfRun, distance, lineRun = False):
    ser.close()
    ser.open()
    ser.write(b"\r\n\r\n")
    time.sleep(2)     # Wait for Printrbot to initialize
    ser.flushInput()  # Flush startup text in serial inp

    string = "\r\n\r\n"
    ser.write(str.encode(string))
    time.sleep(2)     # Wait for Printrbot to initialize
    ser.flushInput()  # Flush startup text in serial inp
    ### start up 
    string = 'M18 ; \n'
    ser.write(str.encode(string))
    time.sleep(2)     # in units of seconds 

    ### set absolute positioning ###
    string = 'G90 ; \n'
    ser.write(str.encode(string))
    time.sleep(2)

    ### set units to mm ###
    string = 'G21 ; \n' 
    ser.write(str.encode(string))
    time.sleep(2)
    
    feedrate = getFeedRate(timeOfRun,distance)
    print("Feedrate: " + str(feedrate))
    print("Distance: " + str(distance) + " mm")
    print("Time Calculated (From Feedrate): " + str(distance/(feedrate/60)))

    i = 0
    points.append(points[0]) # Closes the loop

    User_Input = 1
    while (User_Input == 1):
        for j, point in enumerate(points):

            if i ==0:
                time.sleep(2)
                string = 'G17 X{} Y{} F{} ; \n'.format(centerx,centery,5000)
                ser.write(str.encode(string))
                string = 'G17 X{} Y{} F{} ; \n'.format(centerx+point[0],centery+point[1],5000)
                ser.write(str.encode(string))
                ser.write(b"M400\n")
                time.sleep(0.1)

                wait_for_text('go')
                i += 1
            else:
            # Loop through points on the shape with center at (x=198,y=215)
            # Note that this center point is setup specific and was found by trial-and-error
            # This center point places the center of the shape at the center of the scanner when
            # the 2D stage is mounted onto the top of the scanner 
                #print(j)
                
                if (j < len(points) - 10 and (points[j-1][0] - point[0] > 20 or  points[j-1][0] - point[0] < -20)):
                    #ser.write(b'$$\n')
                    #time.sleep(0.1)

                    # Read the respons 
                    #response = ser.readlines()
                    #print("Response: " + str(response))
                    #print("\n" + "yes" + "\n \n \n \n \n \n")
                    feedrateSPEED = 30000
                    string = 'G1 X{} Y{} F{} ; \n'.format(centerx+point[0],centery+point[1],feedrateSPEED)

                else:
                    string = 'G1 X{} Y{} F{} ; \n'.format(centerx+point[0],centery+point[1],feedrate)

                #print(string)
                ser.write(str.encode(string))

                if(i%20 == 0):
                    time.sleep(((timeOfRun)/len(points))*18)
                i+=1
        if(lineRun):
            for i in range(0, int(timeOfRun/100)):
                time.sleep(99)
                ser.write(b"G4 P0\n")
                
        string = 'G1 X{} Y{} F{} ; \n'.format(centerx+points[0][0],centery+points[0][1],5000)
        ser.write(str.encode(string))
        

            
        user_input_str = '2'
        while (user_input_str != '1' and user_input_str != '0'):
            user_input_str = input("Run Again (1) or Stop (0): ")
            try:
                User_Input = int(user_input_str)
                if (User_Input != 1 and User_Input != 0):
                    print("Enter 1 or 0")
            except ValueError:
                print("Enter an integer")
        i = 1
                
        

    wait_for_text('go')
    string = 'G17 X{} Y{} F5000 ; \n'.format(0,0)
    ser.write(str.encode(string))
    time.sleep(2)
       
def main():
    # Get's run type from user
    user_type_run = "K"
    while (user_type_run != "L" and user_type_run != "P" and user_type_run != "C" and user_type_run != "S" and user_type_run != "LV" and user_type_run != "LH"):
        user_type_run = input("Type of run Smooth (S) or Pausing (P) or Circular (C) or Vertical Line (LV) or Horizontal Line (LH): ")
        if (user_type_run != "LV" and user_type_run != "P" and user_type_run != "C" and user_type_run != "S" and user_type_run != "LH"):
            print("Please enter either L or P or S or C")
    
    # Initializes useful variables
    radius = 0
    line_len = -1
    
    # Gets circle radius/line length from user depending on type of run
    if (user_type_run == "C"):
        while (radius <= 0):
            user_rad_str = input("Enter a radius: ")
            try:
                radius = int(user_rad_str)
            except ValueError:
                print("Please enter an integer")
    elif (user_type_run == "LH" or user_type_run == "LV"):
        while (line_len < 0):
            line_len_str = input("Enter length (mm): ")
            try: 
                line_len = int(line_len_str)
            except ValueError:
                print("Please enter an integer")
    
    # Trace's image from image file and returns total distance that will be travelled by the motor
    if (user_type_run == "P" or user_type_run == "S"):
        print('Pass')
        pass
        #distance = get_source_loc.main()
        #source_loc = read_loc()
    
    # Get's time of run from user
    user_time_num = 0
    while (user_time_num <= 10):
            user_time_num_str = input("Time of run (seconds): ")
            try:
                user_time_num = int(user_time_num_str)
                if (user_time_num <= 10):
                    print("Please enter a time greater than 10 seconds")
            except ValueError:
                print("Please enter an integer")
    
    '''    
    user_type_run = "K"
    while (user_type_run != "L" and user_type_run != "P" and user_type_run != "C"):
        user_type_run = input("Type of run Linear (L) or Pausing (P) or Circular (C): ")
        if (user_type_run != "L" and user_type_run != "P" and user_type_run != "C"):
            print("Please enter either L or P")
            
    radius = 0
    if (user_type_run == "C"):
        while (radius <= 0):
            user_rad_str = input("Enter a radius: ")
            try:
                radius = int(user_rad_str)
            except ValueError:
                print("Please enter an integer")

    '''
    
    # Run's the function that controls the motor
    if(user_type_run == "P"):
        spatialresolutionpoints1cm = [[0, 10], [-7.071, 7.071], [-10, 0], [-7.071, -7.071], [0, -10], [7.071, -7.071], [10, 0], [7.071, 7.071]]
        spatialresolutionpoints10cm = [[0, 100], [-70.71, 70.71], [-100, 0], [-70.71, -70.71], [0, -100], [70.71, -70.71], [100, 0], [70.71, 70.71]]
        toftest = [[0, 0], [5, 0], [10, 0], [20, 0], [40, 0], [80, 0]]
        discretePoints(toftest,user_time_num)
        #shapeMotion(source_loc,user_time_num)
    elif (user_type_run == "S"):
        linearShapeMotion(source_loc,user_time_num,distance)
    elif(user_type_run == "C"): 
        distance = get_source_loc_cir.main(radius) # Get's circle shape
        source_loc = read_loc()
        circularMotion(source_loc,user_time_num,distance)
    elif(user_type_run == "LH"):
        #source_loc = [[0,-line_len/2],[0,line_len/2]] # Get's two points necessary for line
        source_loc = [[line_len/2,0],[-line_len/2,0]]
        #linearShapeMotion(source_loc,user_time_num,line_len, True)
        start_point, end_point = [-line_len/2,0], [line_len/2,0]
        discreteLinearShapeMotion(start_point, end_point, line_len // 10 + 1, user_time_num, True)
    else:
        source_loc = [[0,line_len/2],[0,-line_len/2]]
        linearShapeMotion(source_loc,user_time_num,line_len, True)
    return

main()
# This is the actual line that executes the above functions to move the motor in the circular motion of a given radius in a given time
# First parameter is the radius (in mm) of the circle and the last parameters is run time/duration (in seconds)
#circulatMotion(100, 1800) # Original parameters
# circulatMotion(0, 5) # Test parameters (setting radius to 0 to align the center of the circle with the center of the scanner)
#circulatMotion(60, 7200) # Run parameters 
