import os
import csv
import random
from moviepy.editor import *

os.chdir('/Users/Dustin/Dropbox/DSCN/Experiments/TRW/nonsocial')

f = open('nonsocial_scenes.csv') # CHANGE TO UNIX LINE ENDINGS!!
csv_f = csv.reader(f)

# r = open('Body Bus Random Order.csv')
# csv_r = csv.reader(r)

scenes = []
ons = []
offs = []
for row in csv_f:
    ons.append(row[0])
    offs.append(row[1])
    scenes.append(row[2])

scenes = map(int, scenes)
 
#r_scene = []
#for row in csv_r:
#	r_scene.append(row[0])

#r_scene = scenes  
#random.seed()
#r_scene = random.sample(scenes,len(scenes))
#r_scene = map(int, r_scene)
r_scene = [6,2,9,3,8,11,4,1,12,7,5,10]

##write a random order to analyze
#rand = zip(r_scene)
#with open('/Users/Dustin/Dropbox/DSCN/Experiments/TRW/nonsocial/nonsocial_rand.csv', 'wb') as f:
#	writer = csv.writer(f)
#	writer.writerows(rand)

clip = VideoFileClip('3M Brand Machine.mp4')

temp = []
for row in r_scene:
    print(row)
    temp.append(clip.subclip(ons[row-1], offs[row-1]))
    
catclip = concatenate_videoclips(temp)
catclip_resize = catclip.resize((1024,576))

# catclip.write_videofile('nonsocial_scrambled.mp4',audio=False,threads=20)
catclip.write_videofile('nonsocial_scrambled.mp4',audio=True,remove_temp=False)

