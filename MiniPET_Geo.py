#  MiniPET Geometry

import numpy as np
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

### True if crystal center desired for reconstruction, False if crystal face desired
crystal_center = True

#definition of constants, distances in mm, angles in radians
pi = 3.14159265
radius = 210 #for MiniPET
crystal_center_distance = 7.5 #15.0 / 2

# gap between 64sized arrays
gap = 0.60 #mm
# pixel length each is 3.000mm crystal with 0.200mm spacing between
pix = 3.200 #mm 

# Keeping conventions from Kyle's TPPT_Geo:
# ModCol: not used here, there is only 1 module per side
# ArrCol: 0-7 per module, 0 corresponds to smallest y
# ArrRow: 0-15, 0 corresponds to smallest z

# Rewrote the following functions from Kyle's to be simplified for MiniPET
def max_arr_col(): # number of columns per module
    return 16

def max_arr_row(): # number of rows per module
    return 8

def num_mod(): # number of modules
    return 2

def middle_column_x():
    if crystal_center:
        return radius + crystal_center_distance
    return radius

def middle_column_y():
    return 0

def middle_column_z():
    return 0

def nth_column_x():
    return middle_column_x()

def nth_column_y(ArrCol):
    # Starting from 0, need half a pixel, 7 more modules to the right, and half the gap to reach the middle of the right column
    zeroth_column_y = middle_column_y() + (pix/2+0.1) + 7*pix + gap/2
    return round(zeroth_column_y - (pix*ArrCol + gap*(ArrCol // 8) + 0.2*(ArrCol//8)),1)

def nth_row_z(ArrRow):
    # Starting from 0, need half a pixel, the 3 more modules above that to reach the middle of the top row
    zeroth_column_z = middle_column_z() + (pix/2) + 3*pix
    return round(zeroth_column_z - pix*ArrRow,1)

df = pd.DataFrame(columns=["GeoChannelID","x", "y", "z"])
for i in range(num_mod()): #2 sides left and right
    for k in range(max_arr_col()): #8 ArrCols per module
        for m in range(max_arr_row()): #16 ArrRows in total
            currID = i*max_arr_col()*max_arr_row()+k*max_arr_row()+m
            if i < num_mod()/2:
                df = pd.concat([df, pd.DataFrame([[currID, nth_column_x(), nth_column_y(k), nth_row_z(m)]], columns=["GeoChannelID","x", "y", "z"])], ignore_index=True)
            else:
                df = pd.concat([df, pd.DataFrame([[currID, -1*nth_column_x(), nth_column_y(k), nth_row_z(m)]], columns=["GeoChannelID","x", "y", "z"])], ignore_index=True)

#df = df.sort_values(by="GeoChannelID").reset_index(drop=True)
#df.to_csv("MiniPet_FLASH_Scanner_map_5_19_25.csv", header=False, index=False)