import os
import csv
import random
from moviepy.editor import *

os.chdir('/Users/Dustin/Desktop')

f = open('BB_edited_times.csv') # CHANGE TO UNIX LINE ENDINGS!!
csv_f = csv.reader(f)

scenes = []
ons = []
offs = []
for row in csv_f:
    ons.append(row[0])
    offs.append(row[1])
    scenes.append(row[2])
   
r_scene = scenes  
r_scene = map(int, r_scene)

clip = VideoFileClip("The Body Bus.mkv")

dur = []
for row in r_scene:
    	temp = clip.subclip(ons[row-1], offs[row-1])
    	dur.append(round(temp.duration,4))
    	print(dur)

dur = zip(dur)
with open('Body_Bus_NEW_durations.csv', 'wb') as f:
	writer = csv.writer(f)
	writer.writerows(dur)