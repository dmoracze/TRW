#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy2 Experiment Builder (v1.83.04),
    on Mon Jul 25 21:57:08 2016
If you publish work using this script please cite the PsychoPy publications:
    Peirce, JW (2007) PsychoPy - Psychophysics software in Python.
        Journal of Neuroscience Methods, 162(1-2), 8-13.
    Peirce, JW (2009) Generating stimuli for neuroscience using PsychoPy.
        Frontiers in Neuroinformatics, 2:10. doi: 10.3389/neuro.11.010.2008
"""

from __future__ import absolute_import, division

import psychopy
psychopy.useVersion('latest')

from psychopy import locale_setup, visual, core, data, event, logging, sound, gui, microphone
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions
import sys  # to get file system encoding

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__)).decode(sys.getfilesystemencoding())
os.chdir(_thisDir)

# Store info about the experiment session
expName = 'TRWcomp'  # from the Builder filename that created this script
expInfo = {u'episode': u'', u'participant': u'', u'condition': u''}
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName, order=['participant','episode','condition'])
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s.comp' %(expInfo['participant'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath=None,
    savePickle=True, saveWideText=True,
    dataFileName=filename+'/'+expInfo['participant']+'.'+expInfo['episode']+expInfo['condition'])
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'/'+expInfo['participant']+'.'+expInfo['episode']+expInfo['condition']+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Start Code - component code to be run before the window creation
wavDirName = filename + '/audio'
if not os.path.isdir(wavDirName):
    os.makedirs(wavDirName)  # to hold .wav files

# Setup the Window
win = visual.Window(
    size=(1280, 800), fullscr=True, screen=0,
    allowGUI=False, allowStencil=False,
    monitor='testMonitor', color=[-1,-1,-1], colorSpace='rgb',
    blendMode='avg', useFBO=True)

# Enable sound input/output:
microphone.switchOn()
Qnum = 1
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# Initialize components for Routine "I"
IClock = core.Clock()
instruction = visual.TextStim(win=win, name='instruction',
    text='default text',
    font='Arial',
    units='norm', pos=[0, 0], height=0.1, wrapWidth=1.25, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=0.0);

# Initialize components for Routine "ISI"
ISIClock = core.Clock()
isi = core.StaticPeriod(win=win, screenHz=expInfo['frameRate'], name='isi')

# Initialize components for Routine "P"
PClock = core.Clock()
screen_shot = visual.ImageStim(
    win=win, name='screen_shot',units='pix', 
    image='sin', mask=None,
    ori=0, pos=[0, 0], size=[768,432],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=0.0)

# Initialize components for Routine "ISI"
ISIClock = core.Clock()
isi = core.StaticPeriod(win=win, screenHz=expInfo['frameRate'], name='isi')

# Initialize components for Routine "Q"
QClock = core.Clock()
current_question = visual.TextStim(win=win, name='current_question',
    text='default text',
    font='Arial',
    units='norm', pos=[0, -.58], height=0.1, wrapWidth=1.6, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=0.0);
a_pic = visual.ImageStim(
    win=win, name='a_pic',units='pix', 
    image='./pic/Atticus.png', mask=None,
    ori=0, pos=[-350,250], size=[200,200],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-2.0)
a_text = visual.TextStim(win=win, name='a_text',
    text='Atticus',
    font='Arial',
    units='norm', pos=[-.55,.33], height=.075, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=-3.0);
b_pic = visual.ImageStim(
    win=win, name='b_pic',units='pix', 
    image='./pic/Battie.png', mask=None,
    ori=0, pos=[0,250], size=[200,200],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-4.0)
b_text = visual.TextStim(win=win, name='b_text',
    text='Battie',
    font='Arial',
    units='norm', pos=[0,.33], height=.075, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=-5.0);
m_pic = visual.ImageStim(
    win=win, name='m_pic',units='pix', 
    image='./pic/Melanie.png', mask=None,
    ori=0, pos=[350,250], size=[200,200],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-6.0)
m_text = visual.TextStim(win=win, name='m_text',
    text='Melanie',
    font='Arial',
    units='norm', pos=[.55,.33], height=.075, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=-7.0);
d_pic = visual.ImageStim(
    win=win, name='d_pic',units='pix', 
    image='./pic/DebraJo.png', mask=None,
    ori=0, pos=[-350,0], size=[200,200],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-8.0)
d_text = visual.TextStim(win=win, name='d_text',
    text='Debra Jo',
    font='Arial',
    units='norm', pos=[-.55,-.3], height=.075, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=-9.0);
r_pic = visual.ImageStim(
    win=win, name='r_pic',units='pix', 
    image='./pic/Rory.png', mask=None,
    ori=0, pos=[0,0], size=[200,200],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-10.0)
r_text = visual.TextStim(win=win, name='r_text',
    text='Rory',
    font='Arial',
    units='norm', pos=[0,-.3], height=.075, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=-11.0);
t_pic = visual.ImageStim(
    win=win, name='t_pic',units='pix', 
    image='./pic/Tamara.png', mask=None,
    ori=0, pos=[350,0], size=[200,200],
    color=[1,1,1], colorSpace='rgb', opacity=1,
    flipHoriz=False, flipVert=False,
    texRes=128, interpolate=True, depth=-12.0)
t_text = visual.TextStim(win=win, name='t_text',
    text='Tamara',
    font='Arial',
    units='norm', pos=[.55,-.3], height=.075, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=-13.0);

# Initialize components for Routine "ISI"
ISIClock = core.Clock()
isi = core.StaticPeriod(win=win, screenHz=expInfo['frameRate'], name='isi')

# Initialize components for Routine "end"
endClock = core.Clock()
end_text = visual.TextStim(win=win, name='end_text',
    text='Finished!',
    font='Arial',
    pos=[0, 0], height=0.1, wrapWidth=None, ori=0, 
    color='white', colorSpace='rgb', opacity=1,
    depth=0.0);

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# set up handler to look after randomisation of conditions etc
instructions = data.TrialHandler(nReps=1, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=data.importConditions('other/instructions.csv'),
    seed=None, name='instructions')
thisExp.addLoop(instructions)  # add the loop to the experiment
thisInstruction = instructions.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisInstruction.rgb)
if thisInstruction != None:
    for paramName in thisInstruction.keys():
        exec(paramName + '= thisInstruction.' + paramName)

for thisInstruction in instructions:
    currentLoop = instructions
    # abbreviate parameter names if possible (e.g. rgb = thisInstruction.rgb)
    if thisInstruction != None:
        for paramName in thisInstruction.keys():
            exec(paramName + '= thisInstruction.' + paramName)
    
    # ------Prepare to start Routine "I"-------
    t = 0
    IClock.reset()  # clock
    frameN = -1
    continueRoutine = True
    # update component parameters for each repeat
    instruction.setText(instr)
    go_on_I = event.BuilderKeyResponse()
    # keep track of which components have finished
    IComponents = [instruction, go_on_I]
    for thisComponent in IComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    # -------Start Routine "I"-------
    while continueRoutine:
        # get current time
        t = IClock.getTime()
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *instruction* updates
        if t >= 0.0 and instruction.status == NOT_STARTED:
            # keep track of start time/frame for later
            instruction.tStart = t
            instruction.frameNStart = frameN  # exact frame index
            instruction.setAutoDraw(True)
        
        # *go_on_I* updates
        if t >= 0.0 and go_on_I.status == NOT_STARTED:
            # keep track of start time/frame for later
            go_on_I.tStart = t
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
    
    # -------Ending Routine "I"-------
    for thisComponent in IComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # check responses
    if go_on_I.keys in ['', [], None]:  # No response was made
        go_on_I.keys=None
    instructions.addData('go_on_I.keys',go_on_I.keys)
    if go_on_I.keys != None:  # we had a response
        instructions.addData('go_on_I.rt', go_on_I.rt)
    # the Routine "I" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
# completed 1 repeats of 'instructions'


# ------Prepare to start Routine "ISI"-------
t = 0
ISIClock.reset()  # clock
frameN = -1
continueRoutine = True
routineTimer.add(0.500000)
# update component parameters for each repeat
# keep track of which components have finished
ISIComponents = [isi]
for thisComponent in ISIComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

# -------Start Routine "ISI"-------
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = ISIClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    # *isi* period
    if t >= 0.0 and isi.status == NOT_STARTED:
        # keep track of start time/frame for later
        isi.tStart = t
        isi.frameNStart = frameN  # exact frame index
        isi.start(.5)
    elif isi.status == STARTED:  # one frame should pass before updating params and completing
        isi.complete()  # finish the static period
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in ISIComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "ISI"-------
for thisComponent in ISIComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

# set up handler to look after randomisation of conditions etc
pictures = data.TrialHandler(nReps=1, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=data.importConditions('./pic_order/'+expInfo['episode']+expInfo['condition']+'_Pics.csv'),
    seed=None, name='pictures')
thisExp.addLoop(pictures)  # add the loop to the experiment
thisPicture = pictures.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisPicture.rgb)
if thisPicture != None:
    for paramName in thisPicture.keys():
        exec(paramName + '= thisPicture.' + paramName)

for thisPicture in pictures:
    currentLoop = pictures
    # abbreviate parameter names if possible (e.g. rgb = thisPicture.rgb)
    if thisPicture != None:
        for paramName in thisPicture.keys():
            exec(paramName + '= thisPicture.' + paramName)
    
    # ------Prepare to start Routine "P"-------
    t = 0
    PClock.reset()  # clock
    frameN = -1
    continueRoutine = True
    routineTimer.add(3.000000)
    # update component parameters for each repeat
    screen_shot.setImage('./screenshots/'+pic+'.png')
    key_resp_2 = event.BuilderKeyResponse()
    # keep track of which components have finished
    PComponents = [screen_shot, key_resp_2]
    for thisComponent in PComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    # -------Start Routine "P"-------
    while continueRoutine and routineTimer.getTime() > 0:
        # get current time
        t = PClock.getTime()
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *screen_shot* updates
        if t >= 0.0 and screen_shot.status == NOT_STARTED:
            # keep track of start time/frame for later
            screen_shot.tStart = t
            screen_shot.frameNStart = frameN  # exact frame index
            screen_shot.setAutoDraw(True)
        frameRemains = 0.0 + 3- win.monitorFramePeriod * 0.75  # most of one frame period left
        if screen_shot.status == STARTED and t >= frameRemains:
            screen_shot.setAutoDraw(False)
        
        # *key_resp_2* updates
        if t >= 0.0 and key_resp_2.status == NOT_STARTED:
            # keep track of start time/frame for later
            key_resp_2.tStart = t
            key_resp_2.frameNStart = frameN  # exact frame index
            key_resp_2.status = STARTED
            # keyboard checking is just starting
            win.callOnFlip(key_resp_2.clock.reset)  # t=0 on next screen flip
            event.clearEvents(eventType='keyboard')
        frameRemains = 0.0 + 3- win.monitorFramePeriod * 0.75  # most of one frame period left
        if key_resp_2.status == STARTED and t >= frameRemains:
            key_resp_2.status = STOPPED
        if key_resp_2.status == STARTED:
            theseKeys = event.getKeys(keyList=['g'])
            
            # check for quit:
            if "escape" in theseKeys:
                endExpNow = True
            if len(theseKeys) > 0:  # at least one key was pressed
                key_resp_2.keys = theseKeys[-1]  # just the last key pressed
                key_resp_2.rt = key_resp_2.clock.getTime()
                # a response ends the routine
                continueRoutine = False
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in PComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=["escape"]):
            core.quit()
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "P"-------
    for thisComponent in PComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # check responses
    if key_resp_2.keys in ['', [], None]:  # No response was made
        key_resp_2.keys=None
    pictures.addData('key_resp_2.keys',key_resp_2.keys)
    if key_resp_2.keys != None:  # we had a response
        pictures.addData('key_resp_2.rt', key_resp_2.rt)
# completed 1 repeats of 'pictures'


# ------Prepare to start Routine "ISI"-------
t = 0
ISIClock.reset()  # clock
frameN = -1
continueRoutine = True
routineTimer.add(0.500000)
# update component parameters for each repeat
# keep track of which components have finished
ISIComponents = [isi]
for thisComponent in ISIComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

# -------Start Routine "ISI"-------
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = ISIClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    # *isi* period
    if t >= 0.0 and isi.status == NOT_STARTED:
        # keep track of start time/frame for later
        isi.tStart = t
        isi.frameNStart = frameN  # exact frame index
        isi.start(.5)
    elif isi.status == STARTED:  # one frame should pass before updating params and completing
        isi.complete()  # finish the static period
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in ISIComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "ISI"-------
for thisComponent in ISIComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

# set up handler to look after randomisation of conditions etc
questions = data.TrialHandler(nReps=1, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=data.importConditions('./questions/'+expInfo['episode']+expInfo['condition']+'.csv'),
    seed=None, name='questions')
thisExp.addLoop(questions)  # add the loop to the experiment
thisQuestion = questions.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisQuestion.rgb)
if thisQuestion != None:
    for paramName in thisQuestion.keys():
        exec(paramName + '= thisQuestion.' + paramName)

for thisQuestion in questions:
    currentLoop = questions
    # abbreviate parameter names if possible (e.g. rgb = thisQuestion.rgb)
    if thisQuestion != None:
        for paramName in thisQuestion.keys():
            exec(paramName + '= thisQuestion.' + paramName)
    
    # ------Prepare to start Routine "Q"-------
    t = 0
    QClock.reset()  # clock
    frameN = -1
    continueRoutine = True
    # update component parameters for each repeat
    current_question.setText(question)
    go_on_Q = event.BuilderKeyResponse()
    wavName = wavDirName+'/'+expInfo['participant']+'.'+expInfo['episode']+expInfo['condition']+'_'+str(Qnum)+'.wav'
    record = microphone.AdvAudioCapture(name='record', filename=wavName, stereo=False)
    # keep track of which components have finished
    QComponents = [current_question, go_on_Q, a_pic, a_text, b_pic, b_text, m_pic, m_text, d_pic, d_text, r_pic, r_text, t_pic, t_text, record]
    for thisComponent in QComponents:
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    
    # -------Start Routine "Q"-------
    while continueRoutine:
        # get current time
        t = QClock.getTime()
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *current_question* updates
        if t >= 0.0 and current_question.status == NOT_STARTED:
            # keep track of start time/frame for later
            current_question.tStart = t
            current_question.frameNStart = frameN  # exact frame index
            current_question.setAutoDraw(True)
        
        # *go_on_Q* updates
        if t >= 0.0 and go_on_Q.status == NOT_STARTED:
            # keep track of start time/frame for later
            go_on_Q.tStart = t
            go_on_Q.frameNStart = frameN  # exact frame index
            go_on_Q.status = STARTED
            # keyboard checking is just starting
            win.callOnFlip(go_on_Q.clock.reset)  # t=0 on next screen flip
            event.clearEvents(eventType='keyboard')
        if go_on_Q.status == STARTED:
            theseKeys = event.getKeys(keyList=['space'])
            
            # check for quit:
            if "escape" in theseKeys:
                endExpNow = True
            if len(theseKeys) > 0:  # at least one key was pressed
                go_on_Q.keys = theseKeys[-1]  # just the last key pressed
                go_on_Q.rt = go_on_Q.clock.getTime()
                # a response ends the routine
                continueRoutine = False
        
        # *a_pic* updates
        if t >= 0.0 and a_pic.status == NOT_STARTED:
            # keep track of start time/frame for later
            a_pic.tStart = t
            a_pic.frameNStart = frameN  # exact frame index
            a_pic.setAutoDraw(True)
        
        # *a_text* updates
        if t >= 0.0 and a_text.status == NOT_STARTED:
            # keep track of start time/frame for later
            a_text.tStart = t
            a_text.frameNStart = frameN  # exact frame index
            a_text.setAutoDraw(True)
        
        # *b_pic* updates
        if t >= 0.0 and b_pic.status == NOT_STARTED:
            # keep track of start time/frame for later
            b_pic.tStart = t
            b_pic.frameNStart = frameN  # exact frame index
            b_pic.setAutoDraw(True)
        
        # *b_text* updates
        if t >= 0.0 and b_text.status == NOT_STARTED:
            # keep track of start time/frame for later
            b_text.tStart = t
            b_text.frameNStart = frameN  # exact frame index
            b_text.setAutoDraw(True)
        
        # *m_pic* updates
        if t >= 0.0 and m_pic.status == NOT_STARTED:
            # keep track of start time/frame for later
            m_pic.tStart = t
            m_pic.frameNStart = frameN  # exact frame index
            m_pic.setAutoDraw(True)
        
        # *m_text* updates
        if t >= 0.0 and m_text.status == NOT_STARTED:
            # keep track of start time/frame for later
            m_text.tStart = t
            m_text.frameNStart = frameN  # exact frame index
            m_text.setAutoDraw(True)
        
        # *d_pic* updates
        if t >= 0.0 and d_pic.status == NOT_STARTED:
            # keep track of start time/frame for later
            d_pic.tStart = t
            d_pic.frameNStart = frameN  # exact frame index
            d_pic.setAutoDraw(True)
        
        # *d_text* updates
        if t >= 0.0 and d_text.status == NOT_STARTED:
            # keep track of start time/frame for later
            d_text.tStart = t
            d_text.frameNStart = frameN  # exact frame index
            d_text.setAutoDraw(True)
        
        # *r_pic* updates
        if t >= 0.0 and r_pic.status == NOT_STARTED:
            # keep track of start time/frame for later
            r_pic.tStart = t
            r_pic.frameNStart = frameN  # exact frame index
            r_pic.setAutoDraw(True)
        
        # *r_text* updates
        if t >= 0.0 and r_text.status == NOT_STARTED:
            # keep track of start time/frame for later
            r_text.tStart = t
            r_text.frameNStart = frameN  # exact frame index
            r_text.setAutoDraw(True)
        
        # *t_pic* updates
        if t >= 0.0 and t_pic.status == NOT_STARTED:
            # keep track of start time/frame for later
            t_pic.tStart = t
            t_pic.frameNStart = frameN  # exact frame index
            t_pic.setAutoDraw(True)
        
        # *t_text* updates
        if t >= 0.0 and t_text.status == NOT_STARTED:
            # keep track of start time/frame for later
            t_text.tStart = t
            t_text.frameNStart = frameN  # exact frame index
            t_text.setAutoDraw(True)
        
        # *record* updates
        if t >= 0.0 and record.status == NOT_STARTED:
            # keep track of start time/frame for later
            record.tStart = t
            record.frameNStart = frameN  # exact frame index
            record.status = STARTED
            record.record(sec=100, block=False)  # start the recording thread
        
        if record.status == STARTED and not record.recorder.running:
            record.status = FINISHED
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in QComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # check for quit (the Esc key)
        if endExpNow or event.getKeys(keyList=["escape"]):
            core.quit()
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "Q"-------
    for thisComponent in QComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # check responses
    if go_on_Q.keys in ['', [], None]:  # No response was made
        go_on_Q.keys=None
    questions.addData('go_on_Q.keys',go_on_Q.keys)
    if go_on_Q.keys != None:  # we had a response
        questions.addData('go_on_Q.rt', go_on_Q.rt)
    # record stop & responses
    record.stop()  # sometimes helpful
    Qnum = int(Qnum) + 1
    if not record.savedFile:
        record.savedFile = None
    # store data for questions (TrialHandler)
    questions.addData('record.filename', record.savedFile)
    # the Routine "Q" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
# completed 1 repeats of 'questions'


# ------Prepare to start Routine "ISI"-------
t = 0
ISIClock.reset()  # clock
frameN = -1
continueRoutine = True
routineTimer.add(0.500000)
# update component parameters for each repeat
# keep track of which components have finished
ISIComponents = [isi]
for thisComponent in ISIComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

# -------Start Routine "ISI"-------
while continueRoutine and routineTimer.getTime() > 0:
    # get current time
    t = ISIClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    # *isi* period
    if t >= 0.0 and isi.status == NOT_STARTED:
        # keep track of start time/frame for later
        isi.tStart = t
        isi.frameNStart = frameN  # exact frame index
        isi.start(.5)
    elif isi.status == STARTED:  # one frame should pass before updating params and completing
        isi.complete()  # finish the static period
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in ISIComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "ISI"-------
for thisComponent in ISIComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

# ------Prepare to start Routine "end"-------
t = 0
endClock.reset()  # clock
frameN = -1
continueRoutine = True
# update component parameters for each repeat
go_on_end = event.BuilderKeyResponse()
# keep track of which components have finished
endComponents = [end_text, go_on_end]
for thisComponent in endComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED

# -------Start Routine "end"-------
while continueRoutine:
    # get current time
    t = endClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *end_text* updates
    if t >= 0.0 and end_text.status == NOT_STARTED:
        # keep track of start time/frame for later
        end_text.tStart = t
        end_text.frameNStart = frameN  # exact frame index
        end_text.setAutoDraw(True)
    
    # *go_on_end* updates
    if t >= 0.0 and go_on_end.status == NOT_STARTED:
        # keep track of start time/frame for later
        go_on_end.tStart = t
        go_on_end.frameNStart = frameN  # exact frame index
        go_on_end.status = STARTED
        # keyboard checking is just starting
        win.callOnFlip(go_on_end.clock.reset)  # t=0 on next screen flip
        event.clearEvents(eventType='keyboard')
    if go_on_end.status == STARTED:
        theseKeys = event.getKeys(keyList=['space'])
        
        # check for quit:
        if "escape" in theseKeys:
            endExpNow = True
        if len(theseKeys) > 0:  # at least one key was pressed
            go_on_end.keys = theseKeys[-1]  # just the last key pressed
            go_on_end.rt = go_on_end.clock.getTime()
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

# -------Ending Routine "end"-------
for thisComponent in endComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if go_on_end.keys in ['', [], None]:  # No response was made
    go_on_end.keys=None
thisExp.addData('go_on_end.keys',go_on_end.keys)
if go_on_end.keys != None:  # we had a response
    thisExp.addData('go_on_end.rt', go_on_end.rt)
thisExp.nextEntry()
# the Routine "end" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()
# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'/'+expInfo['participant']+'.'+expInfo['episode']+expInfo['condition']+'.csv')
thisExp.saveAsPickle(filename+'/'+expInfo['participant']+'.'+expInfo['episode']+expInfo['condition'])
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()