###################################################################
##                                                      
## Reads mesh file and converts crimping surface mesh in Cartesian 
## coordinates into a specified cylindrical coordinates 
## Writes steps definetion to bring the non-circular cross section 
## of the crimping surface into a circular through multiple steps
## 
## Reads material data specifically for the current size
##
## 
##
## Author: DMa
## Date: Aug 22nd 2017
##
## 08/23/2017: Fixed a bug on titleRow definition. The list inside 
## a list definition was in the loop. This put a new outside list
## in every iteration which is wrong
###################################################################


import time
import os
import re
import numpy as np
from math import *

os.chdir('P')
filePath = 'P'


# Dimeter when device becomes cylinder (mm)
targetD=10
targetDE = targetD / 25.4

shellDiaArray=[1, 2, 3] # inch for each stent size
stentD=[1, 2, 3] #mm for each stent size
funnel_tip_dia = [1, 2, 3] # inch for each stent size

UMAT_parameter_file = 'materialFileName.txt'

sleeve_thickness = 0.1+0.1 # inch, secondary + primary for cath-load
distalShaft=[0.1, 4, -1] # [diameter, shaft end position in positive, shaft end position in negative]

# number of steps between targetD and minD
step_btw = 5

simulationTimes= [1.0, 0.2, 0.1]
massScale = [5E-7, 3E-7] 
stepDescription='Compress to '

# stent size in mm

# two points defining axial direction (unit vector) of local coordinate sys 
c_a = [0.,0.,0.]
c_b = [0.,0.,1.] 


for n in range(len(stentD)):
  # to hold UMAT parameters
  para_UMAT=[]

  stentSize= str(stentD[n])
  thisFileName= 'thisFileName'
  AbaqusInputDeckName='InputFileName.inp'
  meshFileName='Mesh.inp'

# build the regex as a string
  paraStart_regex = "^\*" + re.escape(stentSize) 
  paraStart = re.compile(paraStart_regex)
    
  rowCount = 0
  with open(UMAT_parameter_file, 'r') as f:
    for line in f:
        if paraStart.match(line)!=None: 
           rowCount +=1
        elif rowCount > 0:
           a=[x.strip() for x in line.split(',')]
           para_UMAT.append(a)
           rowCount +=1
           if rowCount == 3:
              break
##
##  
  
  shellDiameter0=shellDiaArray[n]
  minDiameter = (funnel_tip_dia[n] - sleeve_thickness * 2)
  
  Dc = [shellDiameter0, targetDE]
  Ddisp_Inc = (minDiameter - targetDE) / step_btw
  
  for i in range(step_btw):
    Dc.append(Dc[-1] + Ddisp_Inc)
  
  print stentSize + 'mm stent diameter range (inch)'
  for i in Dc:
    print i
  
  shellDisp=[]
  print stentSize + 'mm sten displacements (inch):'
  for i in xrange(1, len(Dc)):
    shellDisp.append((Dc[i] - Dc[i-1]) / 2)
    print shellDisp[-1]
  
  #############################################################
  ###            Read mesh file to match device nodes        ###
  #############################################################
  
  nodeSetStart=re.compile("^\*NODE, NSET=GLOBAL") # re expression to match keyword
  nodeSetEnd=re.compile("^\*") # re expression to match end of global nset
  nodedevice=re.compile("^\*NSET, NSET=Gdevice, GENERATE")
  
  def process(line): # for each line, extract node# and coordinates and create a map entry 
    a= [ x.strip() for x in line.split(',')]
    nodeCoordinates[str(int(a[0]))] = [float(x) for x in a[1:]] # strip() removes white space
  
  def processdevice(line):
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
      if nodedevice.match(line)!=None:
        flagF=1
        continue
      if flagF==1 and nodeSetEnd.match(line)==None:
        deviceNodeCoord=processdevice(line)
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
  
  for i in range(deviceNodeCoord[0], deviceNodeCoord[1]+1):
    a = np.array(nodeCoordinates[str(i)])
    c = R.dot(a) + o_v
    rho = sqrt(c[0]**2 + c[1]**2)  # convert DOF 1 of cylindrical sys
    r_disp[str(i)] = targetDE / 2. - rho  # total radial disp
  
  #############################################################
  ###         Write radial displacement                    ####
  #############################################################
  
  with open(AbaqusInputDeckName, 'wb') as file:
    file.write('** ***************************************************************************\n')
    file.write('**   Written by : ' + filePath + thisFileName + '\n')
    file.write('**   Author     : DMA\n')
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
    for row in para_UMAT:
      for value in row:
         if value != row[-1]:
            file.write(value + ','),
         else: # avoid write ',' after the last entry
            file.write(value)
      file.write('\n') 
    file.write('*Depvar\n')
    file.write('24,\n')
    file.write('*initial conditions, type=temperature\n')
    file.write('GNiTi, 37.0\n')
    file.write('*****************************************************************\n')
    file.write('****************************************\n')
    file.write('*SURFACE, NAME=distal_shaft, TYPE=REVOLUTION                                                          \n')
    file.write('0.,0.,0.,0.,0.,1.                                                                                     \n')
    file.write('START, ' + str(distalShaft[0]/2) + ',' + str(distalShaft[1]) + '                                      \n')
    file.write('LINE, ' + str(distalShaft[0]/2) + ',' + str(distalShaft[2]) + '                                        \n')
    file.write('*RIGID BODY, ANALYTICAL SURFACE=distal_shaft, REF NODE=10000000                                       \n')
    file.write('*node, nset=refnode                                                                                   \n')
    file.write('10000000,\n')
    file.write('*************************************************************************\n')
    file.write('**Definition of the teflon/elastomer covering\n')
    file.write('*Membrane Section, elset=Gmembrane, material=Mptfe\n')
    file.write('1,\n')
    file.write('*Material, name=Mptfe\n')
    file.write('*hyperelastic\n')
    file.write('1., 0.,\n')
    file.write('*Density\n')
    file.write('1,\n')
    file.write('*****************************************************************\n')
    file.write('*Surface Section, elset=Gshell, density=1.\n')
    file.write('*Surface Section, elset=Gdevice, density=1.\n')
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
    file.write('GNiTi,Gmembrane\n')
    file.write('*nset, nset=gnall\n')
    file.write('GNiTi,Gmembrane\n')
    file.write('**********************************************************************************************\n')
    file.write('**********************************************************************************************\n')
    
    # put some assertion of length of stepDescription ==  length of shellDisp
    for i in range(len(shellDisp)):    
      if i == 0 :  # first step, use 0.4s, define every thing 
        stepNum=i+1
        file.write('** \n')
        stepTitle = stepDescription + str(Dc[stepNum]*25.4) + ' mm'
        file.write('** STEP: Step-' + str(stepNum) + ' ' + stepTitle + '\n')
        file.write('** \n')
        file.write('**\n')
        file.write('*Step,NLGEOM=YES, name=STEP-' + str(stepNum) + '\n')
        file.write(stepTitle + '\n')
        file.write('*Dynamic, Explicit\n')
        file.write(', ' + str(simulationTimes[0]) + '\n')
        file.write('*Bulk Viscosity\n')
        file.write('0.06, 1.2\n')
        file.write('**\n')
        file.write('*variable mass scaling, type=below min, dt=' + str(massScale[0]) + ', freq=10000' + '\n')
        #file.write('*Fixed Mass Scaling, TYPE=below min, dt=' + str(massScale[0]) + ', elset=gniti\n')
        #file.write('*Fixed Mass Scaling, TYPE=below min, dt=' + str(massScale[0]) + ', elset=gmembrane\n')
        file.write('** \n')
        file.write('** BOUNDARY CONDITIONS\n')
        split_num = str(simulationTimes[0]).split('.') # split time into integer and decimal parts
        
        file.write('** Name: Crush_Cylinder Type: Velocity/Angular velocity\n')
        file.write('*Boundary, amplitude=UnitMotion_p' + split_num[1] + ', type=VELOCITY\n')    
        file.write('refnode, 1, 6\n')
        #file.write('Gshell, 1, 1, 0.2\n')
        file.write('Gshell, 1, 3\n')
        #file.write('Gdevice, 1, 1, ' + str(0.1) + '\n')
        file.write('Gdevice, 2, 3\n')
        
        for j in range(deviceNodeCoord[0], deviceNodeCoord[1]+1):
          a = r_disp[str(j)] 
          deviceBC = str(j) + ', ' + '1, 1, ' + str(a) + '\n'
          file.write(deviceBC)
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
        file.write('device_surface, shell_surface\n')
        file.write('device_surface, device_surface\n')      
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
        file.write('*Energy Output, elset=GROI\n')
        file.write('ALLIE,ALLKE\n')
        file.write('**\n')
        file.write('*End Step\n')
        file.write('**********************************************************************************************\n')
        file.write('**********************************************************************************************\n')
      else: 
        stepNum=i+1
        file.write('** \n')
        stepTitle = stepDescription + str(Dc[stepNum]*25.4) + ' mm'
        file.write('** STEP: Step-' + str(stepNum) + ' ' + stepTitle + '\n')
        file.write('** \n')
        file.write('**\n')
        file.write('*Step,NLGEOM=YES, name=STEP-' + str(stepNum) + '\n')
        file.write(stepTitle + '\n')
        file.write('*Dynamic, Explicit\n')
        file.write(', ' + str(simulationTimes[2]) + '\n')
        file.write('*Bulk Viscosity\n')
        file.write('0.06, 1.2\n')
        file.write('**\n')
        file.write('*variable mass scaling, type=below min, dt=' + str(massScale[1]) + ', freq=100000' + '\n')
        #file.write('*Fixed Mass Scaling, TYPE=below min, dt=' + str(massScale[0]) + ', elset=gniti\n')
        #file.write('*Fixed Mass Scaling, TYPE=below min, dt=' + str(massScale[0]) + ', elset=gmembrane\n')
        file.write('** \n')
        file.write('** BOUNDARY CONDITIONS\n')
        split_num = str(simulationTimes[2]).split('.') # split time into integer and decimal parts
        file.write('*Boundary, op=NEW, amplitude=UnitMotion_p' + split_num[1] + ', type=VELOCITY\n') # op=NEW is needed because device as a whole is given a new BC  
        file.write('refnode, 1, 6\n') # *boundary, op=new, removes all BCs from previous step
        file.write('Gdevice, 1, 1, ' + str(shellDisp[i]) + '\n')
        file.write('Gdevice, 2, 3\n')
        file.write('Gaxial_constraint,3,3\n')
        file.write('Gtheta_constraint,2,2\n')
        file.write('GShell,1,3\n')
        #file.write('Gshell, 1, 1, ' + str(shellDisp[i]) + '\n')
        file.write('**\n')
        file.write('*End Step\n')
        file.write('**********************************************************************************************\n')
        file.write('**********************************************************************************************\n')
