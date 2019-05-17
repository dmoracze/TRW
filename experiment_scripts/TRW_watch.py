#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy2 Experiment Builder (v1.83.04), June 18, 2016, at 08:34
If you publish work using this script please cite the relevant PsychoPy publications
  Peirce, JW (2007) PsychoPy - Psychophysics software in Python. Journal of Neuroscience Methods, 162(1-2), 8-13.
  Peirce, JW (2009) Generating stimuli for neuroscience using PsychoPy. Frontiers in Neuroinformatics, 2:10. doi: 10.3389/neuro.11.010.2008
"""

from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import locale_setup, visual, core, data, event, logging, sound, gui
from psychopy.constants import *  # things like STARTED, FINISHED
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import sin, cos, tan, log, log10, pi, average, sqrt, std, deg2rad, rad2deg, linspace, asarray
from numpy.random import random, randint, normal, shuffle, np
from jkpsycho import *
import os  # handy system and path functions
import sys # to get file system encoding

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__)).decode(sys.getfilesystemencoding())
os.chdir(_thisDir)

scanner_coms = ScannerComs(port=3, timeout=0.001, baudrate=19200, keyboard=True)

# Store info about the experiment session
expName = 'TRW_run_episodes'  # from the Builder filename that created this script
expInfo = {u'ID': u'', u'run': u'', u'episode': u'', u'order': u''}
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName, order = ["ID", "episode", "order", "run"])
if dlg.OK == False: core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' %(expInfo['ID'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath=None,
    savePickle=True, saveWideText=False,
    dataFileName=filename)
#save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Start Code - component code to be run before the window creation

# Setup the Window
win = visual.Window(size=(1440, 900), fullscr=True, screen=0, allowGUI=False, allowStencil=False,
    monitor='testMonitor', color=[-1,-1,-1], colorSpace='rgb',
    blendMode='avg', useFBO=True,
    )
# store frame rate of monitor if we can measure it successfully
expInfo['frameRate']=win.getActualFrameRate()
if expInfo['frameRate']!=None:
    frameDur = 1.0/round(expInfo['frameRate'])
else:
    frameDur = 1.0/60.0 # couldn't get a reliable measure so guess

# Initialize components for Routine "I"
IClock = core.Clock()
instruction = visual.TextStim(win=win, ori=0, name='instruction',
    text='default text',    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=1.8,
    color='white', colorSpace='rgb', opacity=1,
    depth=0.0)

# Initialize components for Routine "M"
MClock = core.Clock()
movie = visual.MovieStim3(win=win, name='movie',units='pix', 
    noAudio = False,
    filename='./vid/'+expInfo['episode']+expInfo['order']+'_'+expInfo['run']+'.mp4',
    ori=0, pos=[0, 0], opacity=1,
    size=[1152,648],
    depth=0.0,
    )
text_3 = visual.TextStim(win=win, ori=0, name='text_3',
    text=u'+',    font=u'Arial',
    pos=[0, 0], height=0.1, wrapWidth=None,
    color=u'white', colorSpace='rgb', opacity=1,
    depth=-2.0)

# Initialize components for Routine "F"
FClock = core.Clock()
text_2 = visual.TextStim(win=win, ori=0, name='text_2',
    text=u'+',    font=u'Arial',
    pos=[0, 0], height=0.1, wrapWidth=None,
    color=u'white', colorSpace='rgb', opacity=1,
    depth=0.0)

# Initialize components for Routine "end"
endClock = core.Clock()
text = visual.TextStim(win=win, ori=0, name='text',
    text='Finished!',    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=None,
    color='white', colorSpace='rgb', opacity=1,
    depth=0.0)

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# set up handler to look after randomisation of conditions etc
instructions = data.TrialHandler(nReps=1, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=data.importConditions('./other/instructions.csv'),
    seed=None, name='instructions')
thisExp.addLoop(instructions)  # add the loop to the experiment
thisInstruction = instructions.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb=thisInstruction.rgb)
if thisInstruction != None:
    for paramName in thisInstruction.keys():
        exec(paramName + '= thisInstruction.' + paramName)

for thisInstruction in instructions:
    currentLoop = instructions
    # abbreviate parameter names if possible (e.g. rgb = thisInstruction.rgb)
    if thisInstruction != None:
        for paramName in thisInstruction.keys():
            exec(paramName + '= thisInstruction.' + paramName)
    
    #------Prepare to start Routine "I"-------
    t = 0
    IClock.reset()  # clock 
    frameN = -1
    # update component parameters for each repeat
    instruction.setText(instr)
    go_on_I = event.BuilderKeyResponse()  # create an object of type KeyResponse
    go_on_I.status = NOT_STARTED
    # keep track of which components have finished
    IComponents = []
    IComponents.append(instruction)
    IComponents.append(go_on_I)
    for thisComponent in IComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    #-------Start Routine "I"-------
    continueRoutine = True
    while continueRoutine:
        # get current time
        t = IClock.getTime()
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *instruction* updates
        if t >= 0.0 and instruction.status == NOT_STARTED:
            # keep track of start time/frame for later
            instruction.tStart = t  # underestimates by a little under one frame
            instruction.frameNStart = frameN  # exact frame index
            instruction.setAutoDraw(True)
        
        # *go_on_I* updates
        if t >= 0.0 and go_on_I.status == NOT_STARTED:
            # keep track of start time/frame for later
            go_on_I.tStart = t  # underestimates by a little under one frame
            go_on_I.frameNStart = frameN  # exact frame index
            go_on_I.status = STARTED
            # keyboard checking is just starting
            win.callOnFlip(go_on_I.clock.reset)  # t=0 on next screen flip
            event.clearEvents(eventType='keyboard')
        if go_on_I.status == STARTED:
            theseKeys = event.getKeys(keyList=['space'])
            
            # check for quit:
            if "escape" in theseKeys:
                endExpNow = True
            if len(theseKeys) > 0:  # at least one key was pressed
                go_on_I.keys = theseKeys[-1]  # just the last key pressed
                go_on_I.rt = go_on_I.clock.getTime()
                # a response ends the routine
                continueRoutine = False
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in IComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=["escape"]):
            core.quit()
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
            
        if "6" in scanner_coms.messages():
            break
    
    #-------Ending Routine "I"-------
    for thisComponent in IComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # check responses
    if go_on_I.keys in ['', [], None]:  # No response was made
       go_on_I.keys=None
    # store data for instructions (TrialHandler)
    instructions.addData('go_on_I.keys',go_on_I.keys)
    if go_on_I.keys != None:  # we had a response
        instructions.addData('go_on_I.rt', go_on_I.rt)
    # the Routine "I" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
# completed 1 repeats of 'instructions'


#------Prepare to start Routine "M"-------
t = 0
MClock.reset()  # clock 
frameN = -1
# update component parameters for each repeat
go_on_m = event.BuilderKeyResponse()  # create an object of type KeyResponse
go_on_m.status = NOT_STARTED
# keep track of which components have finished
MComponents = []
MComponents.append(movie)
MComponents.append(go_on_m)
MComponents.append(text_3)
for thisComponent in MComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

#-------Start Routine "M"-------
continueRoutine = True
while continueRoutine:
    # get current time
    t = MClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *movie* updates
    if t >= 10 and movie.status == NOT_STARTED:
        # keep track of start time/frame for later
        movie.tStart = t  # underestimates by a little under one frame
        movie.frameNStart = frameN  # exact frame index
        movie.setAutoDraw(True)
    if movie.status == FINISHED:  # force-end the routine
        continueRoutine = False
    
    # *go_on_m* updates
    if t >= 0.0 and go_on_m.status == NOT_STARTED:
        # keep track of start time/frame for later
        go_on_m.tStart = t  # underestimates by a little under one frame
        go_on_m.frameNStart = frameN  # exact frame index
        go_on_m.status = STARTED
        # keyboard checking is just starting
        win.callOnFlip(go_on_m.clock.reset)  # t=0 on next screen flip
        event.clearEvents(eventType='keyboard')
    if go_on_m.status == STARTED:
        theseKeys = event.getKeys(keyList=['g'])
        
        # check for quit:
        if "escape" in theseKeys:
            endExpNow = True
        if len(theseKeys) > 0:  # at least one key was pressed
            go_on_m.keys = theseKeys[-1]  # just the last key pressed
            go_on_m.rt = go_on_m.clock.getTime()
            # a response ends the routine
            continueRoutine = False
    
    # *text_3* updates
    if t >= 0.0 and text_3.status == NOT_STARTED:
        # keep track of start time/frame for later
        text_3.tStart = t  # underestimates by a little under one frame
        text_3.frameNStart = frameN  # exact frame index
        text_3.setAutoDraw(True)
    if text_3.status == STARTED and t >= (0.0 + (10-win.monitorFramePeriod*0.75)): #most of one frame period left
        text_3.setAutoDraw(False)
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in MComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

#-------Ending Routine "M"-------
for thisComponent in MComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if go_on_m.keys in ['', [], None]:  # No response was made
   go_on_m.keys=None
# store data for thisExp (ExperimentHandler)
thisExp.addData('go_on_m.keys',go_on_m.keys)
if go_on_m.keys != None:  # we had a response
    thisExp.addData('go_on_m.rt', go_on_m.rt)
thisExp.nextEntry()
# the Routine "M" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

#------Prepare to start Routine "F"-------
t = 0
FClock.reset()  # clock 
frameN = -1
routineTimer.add(20.000000)
# update component parameters for each repeat
# keep track of which components have finished
FComponents = []
FComponents.append(text_2)
for thisComponent in FComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

#-------Start Routine "F"-------
continueRoutine = True
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = FClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text_2* updates
    if t >= 0.0 and text_2.status == NOT_STARTED:
        # keep track of start time/frame for later
        text_2.tStart = t  # underestimates by a little under one frame
        text_2.frameNStart = frameN  # exact frame index
        text_2.setAutoDraw(True)
    if text_2.status == STARTED and t >= (0.0 + (20-win.monitorFramePeriod*0.75)): #most of one frame period left
        text_2.setAutoDraw(False)
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in FComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

#-------Ending Routine "F"-------
for thisComponent in FComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

#------Prepare to start Routine "end"-------
t = 0
endClock.reset()  # clock 
frameN = -1
# update component parameters for each repeat
go_on_e = event.BuilderKeyResponse()  # create an object of type KeyResponse
go_on_e.status = NOT_STARTED
# keep track of which components have finished
endComponents = []
endComponents.append(text)
endComponents.append(go_on_e)
for thisComponent in endComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

#-------Start Routine "end"-------
continueRoutine = True
while continueRoutine:
    # get current time
    t = endClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text* updates
    if t >= 0.0 and text.status == NOT_STARTED:
        # keep track of start time/frame for later
        text.tStart = t  # underestimates by a little under one frame
        text.frameNStart = frameN  # exact frame index
        text.setAutoDraw(True)
    
    # *go_on_e* updates
    if t >= 0.0 and go_on_e.status == NOT_STARTED:
        # keep track of start time/frame for later
        go_on_e.tStart = t  # underestimates by a little under one frame
        go_on_e.frameNStart = frameN  # exact frame index
        go_on_e.status = STARTED
        # keyboard checking is just starting
        win.callOnFlip(go_on_e.clock.reset)  # t=0 on next screen flip
        event.clearEvents(eventType='keyboard')
    if go_on_e.status == STARTED:
        theseKeys = event.getKeys(keyList=['space'])
        
        # check for quit:
        if "escape" in theseKeys:
            endExpNow = True
        if len(theseKeys) > 0:  # at least one key was pressed
            go_on_e.keys = theseKeys[-1]  # just the last key pressed
            go_on_e.rt = go_on_e.clock.getTime()
            # a response ends the routine
            continueRoutine = False
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in endComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

#-------Ending Routine "end"-------
for thisComponent in endComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if go_on_e.keys in ['', [], None]:  # No response was made
   go_on_e.keys=None
# store data for thisExp (ExperimentHandler)
thisExp.addData('go_on_e.keys',go_on_e.keys)
if go_on_e.keys != None:  # we had a response
    thisExp.addData('go_on_e.rt', go_on_e.rt)
thisExp.nextEntry()
# the Routine "end" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()
# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort() # or data files will save again on exit
win.close()
core.quit()
