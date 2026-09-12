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
import get_source_loc
import get_source_loc_cir
from datetime import datetime

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
                string = 'G17 X{} Y{} F{} ; \n'.format(198,214,5000)
                ser.write(str.encode(string))
                string = 'G17 X{} Y{} F{} ; \n'.format(198+point[0],214+point[1],5000)
                ser.write(str.encode(string))
                wait_for_text('go')
                i += 1
            else:
            # Loops through points on shape
            # Note that this center point is setup specific and was found by trial-and-error
            # This center point places the center of the circle at the center of the scanner when
            # the 2D stage is mounted onto the top of the scanner 
                string = 'G17 X{} Y{} F100 ; \n'.format(198+point[0],214+point[1])
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

def linearShapeMotion(points, timeOfRun, distance):
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
#    points.append(points[0]) # Closes the loop

    User_Input = 1
    while (User_Input == 1):
        for point in points:

            if i ==0:
                time.sleep(2)
                string = 'G17 X{} Y{} F{} ; \n'.format(198,214,5000)
                ser.write(str.encode(string))
                string = 'G17 X{} Y{} F{} ; \n'.format(198+point[0],214+point[1],5000)
                ser.write(str.encode(string))
                ser.write(b"M400\n")
                time.sleep(0.1)

                wait_for_text('go')
                i += 1
            else:
            # Loop through points on the shape with center at (x=198,y=214)
            # Note that this center point is setup specific and was found by trial-and-error
            # This center point places the center of the shape at the center of the scanner when
            # the 2D stage is mounted onto the top of the scanner 

                string = 'G1 X{} Y{} F{} ; \n'.format(198+point[0],214+point[1],feedrate)
                
                ser.write(str.encode(string))

                if(i%20 == 0):
                    time.sleep(((timeOfRun)/len(points))*18)
                i+=1
                
        string = 'G1 X{} Y{} F{} ; \n'.format(198+points[0][0],214+points[0][1],1000)
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


def read_loc():
    data_list = []
    with open(script_dir + "/source_locations", 'r') as file:
        for line in file:
            row = list(map(float, line.strip().split()))
            data_list.append(row)
    return data_list
    
    
def circulatMotion(radius, timeOfRun):
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

    r = radius 
    theta = np.linspace(0,2*np.pi,1000)
    times = getTime(timeOfRun,theta)
    i = 0
    string = 'G17 X{} Y{} F5000 ; \n'.format(198,214)
    ser.write(str.encode(string))
    User_Input = 1
    while (User_Input == 1):
        for point in theta:

            # Loop through points on the circle with center at (x=198,y=214)
            # Note that this center point is setup specific and was found by trial-and-error
            # This center point places the center of the circle at the center of the scanner when
            # the 2D stage is mounted onto the top of the scanner 
            string = 'G17 X{} Y{} F100 ; \n'.format(198+r*np.cos(point),214+r*np.sin(point))
            ser.write(str.encode(string))
            
            if i ==0:
                time.sleep(2)
                wait_for_text('go')
                i += 1
            else:
                time.sleep(times)

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
    
    
    
def main():
    # Get's run type from user
    user_type_run = "K"
    while (user_type_run != "L" and user_type_run != "P" and user_type_run != "C" and user_type_run != "S"):
        user_type_run = input("Type of run Smooth (S) or Pausing (P) or Circular (C) or Line (L): ")
        if (user_type_run != "L" and user_type_run != "P" and user_type_run != "C" and user_type_run != "S"):
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
    elif (user_type_run == "L"):
        while (line_len < 0):
            line_len_str = input("Enter length (mm): ")
            try: 
                line_len = int(line_len_str)
            except ValueError:
                print("Please enter an integer")
    
    # Trace's image from image file and returns total distance that will be travelled by the motor
    if (user_type_run == "P" or user_type_run == "S"):
        distance = get_source_loc.main()
        source_loc = read_loc()
    
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
        shapeMotion(source_loc,user_time_num)
    elif (user_type_run == "S"):
        linearShapeMotion(source_loc,user_time_num,distance)
    elif(user_type_run == "C"): 
        distance = get_source_loc_cir.main(radius) # Get's circle shape
        source_loc = read_loc()
        linearShapeMotion(source_loc,user_time_num,distance)
    else:
        source_loc = [[0,-line_len/2],[0,line_len/2]] # Get's two points necessary for line
        linearShapeMotion(source_loc,user_time_num,line_len)
    return

main()
# This is the actual line that executes the above functions to move the motor in the circular motion of a given radius in a given time
# First parameter is the radius (in mm) of the circle and the last parameters is run time/duration (in seconds)
#circulatMotion(100, 1800) # Original parameters
# circulatMotion(0, 5) # Test parameters (setting radius to 0 to align the center of the circle with the center of the scanner)
#circulatMotion(60, 7200) # Run parameters 
