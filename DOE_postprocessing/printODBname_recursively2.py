# https://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python
import fnmatch
import os
from odbAccess import *
from abaqusConstants import *
import re


def tryint(s):
   try:
       return int(s)
   except:
       return s

def alphanum_key(s):
    return [ tryint(c) for c in re.split('([0-9]+)',s) ]

def sort_naturally(l):
   l.sort(key=alphanum_key)


odbPath='W:\Technical\CTS\PHX HPCC\Launch 02\cts\llargura\DOE_27mm'
#odbPath='C:\AbaqusJobs'
#matches=['C:/AbaqusJobs/53_crush_maxShaftOD_run_v2018.odb']
matches=[]

for root, dirnames, filenames in os.walk(odbPath):
    sort_naturally(dirnames)
    for filename in fnmatch.filter(filenames, '*.odb'):
              matches.append(os.path.join(root, filename))
	 #matches.append(filename)
for i in matches:
    print i
