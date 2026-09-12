
import tkinter as tk
from tkinter.messagebox import askyesno
import platform
import serial
import serial.tools.list_ports
from serial.tools import list_ports
import numpy as np
import os
import time
from datetime import datetime

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

string = 'G17 X{} Y{} F5000 ; \n'.format(0,0)
ser.write(str.encode(string))
