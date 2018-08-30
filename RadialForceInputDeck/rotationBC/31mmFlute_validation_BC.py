import time
import os
import re
import numpy as np
from math import *

# validation: pull the flute to circle in first step and compress as cylinder, for comparison
# with result from a real cylinder

#os.chdir('C:\\Users\\dma\\Documents\\rotationBC')
os.chdir('P:\\ACMDS2_analysis\\Scripts\\RadialForceInputDeck\\rotationBC')

# Dimeter when flute becomes cylinder (mm)
targetD_position=0 # 0-based step number when flute is cylinder

# stent size in mm
stentDia=31

# two points defining axial direction (unit vector) of local coordinate sys 
c_a = [0.,0.,0.]
c_b = [0.,0.,1.] 


#shellDiaArray=[1.2638, 1.3848, 1.5048, 1.6248, 1.8218, 1.9788, 2.1388] # inch
#shellDiametermm=[31,34,37,40,45,49,53] #mm
shellDiaArray=1.2638 # inch
stentDiametermm=31

thisFileName= str(stentDia) + 'mmFlute_BC.py'
AbaqusInputDeckName=str(stentDia) + 'mmFlute_BC_stp' + str(targetD_position+1) + '_validation.inp'
#meshFileName='ACMDS2_' + str(stentDiametermm) + 'mm_crushModel_mesh.inp'
# 1) for reading flute mesh, 2) for print included mesh file name in run file, can be modified after run file is generated
meshFileName = 'ACMDS2_31mm_shorter_coarserStent_mesh.inp' 
#for j in range(len(shellDiametermm)):
shellDiameter0=shellDiaArray
simulationTimes= [0.4, 0.2, 0.1]
#timeStep = 0.4 # step displacement 0.03937 is small
massScale = 8E-7
stepDescription='Compress to '


metricD = shellDiameter0 * 25.4
secondD = floor(metricD)
halfD = floor(metricD / 2.0)
halfD = halfD - halfD % 2 # if odd number plus 1
stepCount = (secondD - halfD) / 2
Dc=[]
Dc.append(metricD)

# range of diameters for output frames
for i in range(int(secondD), int(halfD-1), -2):
  Dc.append(i)

DcE=[]  
print 'Diameter range (mm)'
for i in Dc: 
  DcE.append(i)
  print i
print('**************\n**************\n')

for i in range(len(Dc)):
  Dc[i] = Dc[i]/25.4  

print 'Diameter range (inch)'

for i in Dc: 
  print i  
print('**************\n**************\n')


totalDisp_R = (Dc[targetD_position] - shellDiameter0) /2.  # total radial disp before flute becoming cylindr(negative)
shellDispRatio=[]
shellDispRatio.append( (Dc[0] - shellDiameter0) / 2.0) # ratio of total radial disp for step1

for i in xrange(1, len(Dc), 1):
  shellDispRatio.append((Dc[i] - Dc[i-1]) / 2.0)  # ratio is positive

print('displacement:')
for i in xrange(0, (len(Dc)), 1):
  print shellDispRatio[i],


#############################################################
###            Read mesh file to match flute nodes        ###
#############################################################

nodeSetStart=re.compile("^\*NODE, NSET=GLOBAL") # re expression to match keyword
nodeSetEnd=re.compile("^\*") # re expression to match end of global nset
nodeFlute=re.compile("^\*NSET, NSET=GFlute, GENERATE")

def process(line): # for each line, extract node# and coordinates and create a map entry 
  a= [ x.strip() for x in line.split(',')]
  nodeCoordinates[str(int(a[0]))] = [float(x) for x in a[1:]] # strip() removes white space

def processFlute(line):
  a=[ int(y) for y in [x.strip() for x in line.split(',')]]
  return a

nodeCoordinates={} # map to store all nodes' info
flagG=0
flagF=0
with open(meshFileName) as f:
  for line in f:
    if nodeSetStart.match(line)!=None:
      flagG=1  # turn on flag to process
      continue # get a new line with node 1
    if flagG==1 and nodeSetEnd.match(line)==None:
      process(line) 
    elif flagG==1 and nodeSetEnd.match(line)!=None:
      flagG=0  # nset ends, so turn off flag
    if nodeFlute.match(line)!=None:
      flagF=1
      continue
    if flagF==1 and nodeSetEnd.match(line)==None:
      fluteNodeCoord=processFlute(line)
      falgF=0
      break
#############################################################
###            Calculate rotation matrix                  ###
#############################################################

# this is an implementation of method stated in:
# https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d

# axial unit vector of local coordinate sys
_a = [(y-x) for x,y in zip(c_a,c_b)]

# axial unit vector of global sys
_b = [0.,0.,1.]

# cross product
v=np.cross(_a, _b)
#v=[_a[1] * _b[2] - _a[2] * _b[1],  
#_a[2] * _b[0] - _a[0] * _b[2],
#_a[0] * _b[1] - _a[1] * _b[0]]


# ||v|| (sine of angle)
s = sqrt( v[0]**2 + v[1]**2 + v[2]**2 )

# dot product (cosine of angle)
c =np.dot(_a, _b)
# c=_a[0]*_b[0] + _a[1]*_b[1] + _a[2]*_b[2]

# skew-symmetric cross-product matrix of v
vx = np.array([(0., v[2], -v[1]), (-v[2], 0., v[0]), (v[1], -v[0], 0.)])

# (1-c)/s^2
cof = 1./(1. + c)

# vx^2
vx2 = np.dot(vx,vx) 

# create I
I = np.eye(3)

#
R = I + vx + cof * vx2

#############################################################
###     Rotate and translate global coordinates to local  ###
#############################################################
# to convert global coordinates to local ones, rotate the global coord in opposite angular direction of the angle from global axial # vector to local axial vector. To make these two vectors align each other, also need to translate global coord by the vector from  
# local origin to global origin

o_v = [-x for x in c_a]  # global origin is at (0, 0, 0), so translatoin vector is -c_a

r_disp={}

for i in range(fluteNodeCoord[0], fluteNodeCoord[1]+1):
  a = np.array(nodeCoordinates[str(i)])
  c = R.dot(a) + o_v
  rho = sqrt(c[0]**2 + c[1]**2)  # convert DOF 1 of cylindrical sys
  r_disp[str(i)] = Dc[targetD_position] / 2. - rho  # total radial disp

#############################################################
###         Write radial displacement                    ####
#############################################################

with open(AbaqusInputDeckName, 'wb') as file:
  file.write('** ***************************************************************************\n')
  file.write('**   Written by : ' + thisFileName + '\n')
  file.write('**   Author     : DMA')
  file.write('**   Date       : ' + time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()) + '\n')
  file.write('** ***************************************************************************\n** \n')
  file.write('*include, input=' + meshFileName + '\n')
  file.write('**\n** Load Step 1 -------------------------------------------------------\n')
  file.write('*transform, type=C, nset=global\n')
  file.write('0.,0.,0.,0.,0.,1.\n')
  file.write('*******************************************\n')
  file.write('*tie, name=ptfe_2_stent, position tolerance=0.003, adjust=no, no thickness, type=node to surface\n')
  file.write('cs2, cs1\n')
  file.write('*******************************************\n')
  file.write('** ELEMENT CONTROLS\n**\n') 
  file.write('*Section Controls, name=EC-1, hourglass=ENHANCED\n')
  file.write('1., 1., 1.\n')
    
  # generate amplitudes  
  for i in simulationTimes:
    amp = 2.0/i
    split_num = str(i).split('.') # split time into integer and decimal parts
    file.write('*Amplitude, name=UnitMotion_p' + split_num[1] + ', definition=SMOOTH STEP' + '\n')
    file.write('., 0., ' + str(i/2) + ', ' + str(amp) + ', ' + str(i) + ', 0.\n')
  
  file.write('*Amplitude, name=DSLOAD_AMP, definition=SMOOTH STEP\n')
  file.write('0., 0., ' + str(simulationTimes[0]) + ', 1.\n')
  
  file.write('** \n')
  file.write('**************************************************************************\n')
  file.write('******Solid Section for the Nitinol Wire**********************************\n')
  file.write('*SOLID SECTION, ELSET=GNiTi, controls=EC-1, MATERIAL=abq_super_elastic_n3d\n')
  file.write('*material, name=abq_super_elastic_n3d\n')
  file.write('*Density\n')
  file.write('0.000611979,\n')
  file.write('*user material,constants=14\n')
  file.write('6.21e6,0.35, 3.58e6,0.35,0.0372,206.,77.595e3,84.122e3\n')
  file.write('37.,1142.,49313.,40611,116.4e3,0.0372\n')
  file.write('*Depvar\n')
  file.write('24,\n')
  file.write('*initial conditions, type=temperature\n')
  file.write('GNiTi, 37.0\n')
  file.write('*****************************************************************\n')
  file.write('****************************************\n')
  #file.write('*SURFACE, NAME=distal_shaft, TYPE=REVOLUTION                                                          \n')
  #file.write('0.,0.,0.,0.,0.,1.                                                                                     \n')
  #file.write('START, ' + str(distalShaft[0]/2) + ',' + str(distalShaft[1]) + '                                      \n')
  #file.write('LINE, ' + str(distalShaft[0]/2) + ',' + str(distalShaft[2]) + '                                        \n')
  #file.write('*RIGID BODY, ANALYTICAL SURFACE=distal_shaft, REF NODE=10000000                                       \n')
  #file.write('*node, nset=refnode                                                                                   \n')
  #file.write('10000000,\n')
  #file.write('*************************************************************************\n')
  file.write('**Definition of the teflon/elastomer covering\n')
  file.write('*Membrane Section, elset=Ggraft, material=Mptfe\n')
  file.write('0.0039,\n')
  file.write('*Material, name=Mptfe\n')
  file.write('*hyperelastic\n')
  file.write('2894., 0.,\n')
  file.write('*Density\n')
  file.write('0.00015,\n')
  file.write('*****************************************************************\n')
  file.write('*Surface Section, elset=Gshell, density=1.\n')
  file.write('*Surface Section, elset=Gflute, density=1.\n')
  file.write('*****************************************************************\n')
  file.write('** INTERACTION PROPERTIES\n')
  file.write('** \n')
  file.write('*Surface Interaction, name=cp1\n')
  file.write('*Friction\n')
  file.write('0.01,\n')
  file.write('*Surface Behavior, pressure-overclosure=SCALE FACTOR\n')
  file.write('0.1, ,1.5,\n')
  file.write('****************************************\n')
  file.write('*elset, elset=geall\n')
  file.write('GNiTi,Ggraft\n')
  file.write('*nset, nset=gnall\n')
  file.write('GNiTi,Ggraft\n')
  file.write('**********************************************************************************************\n')
  file.write('**********************************************************************************************\n')
  
  # put some assertion of length of stepDescription ==  length of shellDisp
  for i in xrange(0, len(shellDispRatio), 1):    
    if i == 0 :  # first step, use 0.4s, define every thing 
      stepNum=i+1
      simTime = str(simulationTimes[0])
      file.write('** \n')
      stepTitle = stepDescription + str(DcE[i]) + 'mm'
      file.write('** STEP: Step-' + str(stepNum) + ' ' + stepTitle + '\n')
      file.write('** \n')
      file.write('**\n')
      file.write('*Step,NLGEOM=YES, name=STEP-' + str(stepNum) + '\n')
      file.write(stepTitle + '\n')
      file.write('*Dynamic, Explicit\n')
      file.write(', ' + simTime + '\n')
      file.write('*Bulk Viscosity\n')
      file.write('0.06, 1.2\n')
      file.write('**\n')
      file.write('*Fixed Mass Scaling, TYPE=below min, dt=' + str(massScale) + ', elset=gniti\n')
      file.write('*Fixed Mass Scaling, TYPE=below min, dt=' + str(massScale) + ', elset=ggraft\n')
      file.write('** \n')
      file.write('** BOUNDARY CONDITIONS\n')
      split_num = simTime.split('.') # split time into integer and decimal parts
      
      file.write('** Name: Crush_Cylinder Type: Velocity/Angular velocity\n')
      file.write('*Boundary, amplitude=UnitMotion_p' + split_num[1] + ', type=VELOCITY\n')    
      #file.write('refnode, 1, 6\n')
      file.write('Gshell, 1, 3\n')
      #file.write('Gflute, 1, 1, ' + str(0.1) + '\n')
      file.write('Gflute, 2, 3\n')
      
      for j in range(fluteNodeCoord[0], fluteNodeCoord[1]+1):
        a = r_disp[str(j)] 
        fluteBC = str(j) + ', ' + '1, 1, ' + str(a) + '\n'
        file.write(fluteBC)
      file.write('Gaxial_constraint,3,3\n')
      file.write('Gtheta_constraint,2,2\n')
      file.write('** \n')
      file.write('** \n')
      #file.write('*DSLOAD, amplitude=DSLOAD_AMP\n')
      #file.write('cs3, P, 1.1969\n')
      file.write('** \n')
      file.write('** INTERACTIONS\n')
      file.write('** \n')
      file.write('** Interaction: cp1\n')
      file.write('*Contact, op=NEW\n')
      file.write('*Contact Inclusions, ALL EXTERIOR\n')
      file.write('*Contact Exclusions\n')
      file.write('flute_surface, shell_surface\n')
      file.write('flute_surface, flute_surface\n')
      file.write('*Contact property assignment                \n')
      file.write(' ,  , cp1\n')
      file.write('** \n')
      file.write('** OUTPUT REQUESTS\n')
      file.write('** \n')
      file.write('** \n')
      file.write('** FIELD OUTPUT: F-Output-1\n')
      file.write('** \n')
      file.write('*Output, field, number interval=1\n')
      file.write('*Contact Output\n')
      file.write('cforce\n')
      file.write('*Element Output, elset=Gniti\n')
      file.write('LE,\n')
      file.write('*node output\n')
      file.write('U, RF\n')
      file.write('** \n')
      file.write('**\n')
      file.write('** HISTORY OUTPUT: H-Output-1\n')
      file.write('** \n')
      file.write('*Output, history, time interval=0.02\n')
      file.write('*Energy Output, elset=geall\n')
      file.write('ALLIE,ALLKE\n')
      file.write('*Energy Output, elset=gniti\n')
      file.write('ALLIE,ALLKE\n')
      file.write('**\n')
      file.write('*End Step\n')
      file.write('**********************************************************************************************\n')
      file.write('**********************************************************************************************\n')
    else: 
      stepNum=i+1
      simTime = str(simulationTimes[1])
      file.write('** \n')
      stepTitle = stepDescription + str(DcE[i]) + 'mm'
      file.write('** STEP: Step-' + str(stepNum) + ' ' + stepTitle + '\n')
      file.write('** \n')
      file.write('**\n')
      file.write('*Step,NLGEOM=YES, name=STEP-' + str(stepNum) + '\n')
      file.write(stepTitle + '\n')
      file.write('*Dynamic, Explicit\n')
      file.write(', ' + simTime + '\n')
      file.write('*Bulk Viscosity\n')
      file.write('0.06, 1.2\n')
      file.write('**\n')
      file.write('*Fixed Mass Scaling, TYPE=below min, dt=' + str(massScale) + ', elset=gniti\n')
      file.write('*Fixed Mass Scaling, TYPE=below min, dt=' + str(massScale) + ', elset=ggraft\n')
      file.write('** \n')
      file.write('** BOUNDARY CONDITIONS\n')
      split_num = simTime.split('.') # split time into integer and decimal parts

      # validation has set targetD_position = 0, so this never executes
      if i <= targetD_position:            
        file.write('*Boundary, amplitude=UnitMotion_p' + split_num[1] + ', type=VELOCITY\n')
        for j in range(fluteNodeCoord[0], fluteNodeCoord[1]+1):
          a = shellDispRatio[i] * r_disp[str(j)] 
          fluteBC = str(j) + ', ' + '1, 1, ' + str(a) + '\n'
          file.write(fluteBC)
      else:
	  # op=NEW to remove single node BC in previous to allow flute BC as a whole, but all BCs have to be redefined here.
        file.write('*Boundary, op=NEW, amplitude=UnitMotion_p' + split_num[1] + ', type=VELOCITY\n') 
        a = shellDispRatio[i]
        file.write('GFLUTE, 1, 1, ' + str(a) + '\n')  
        file.write('GFLUTE, 2, 3\n')
        file.write('Gaxial_constraint,3,3\n')
        file.write('Gtheta_constraint,2,2\n')
        file.write('GShell,1,3\n')
      #file.write('Gshell, 1, 1, ' + str(shellDisp[i]) + '\n')
      file.write('**\n')
      file.write('*End Step\n')
      file.write('**********************************************************************************************\n')
      file.write('**********************************************************************************************\n')
