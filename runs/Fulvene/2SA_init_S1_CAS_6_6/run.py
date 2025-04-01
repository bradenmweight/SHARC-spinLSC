#!/usr/bin/python

import random
import os
import shutil
import numpy as np
import subprocess as sp
off = 100
nproc = 1000

memory = 20000 # In words: 1 word = 8 bytes

global NAtoms
NAtoms = 12
NStates = 2

zeroMomentum = True

def getxyz(n,ini = "initconds") :
	n += 1
	fob = [at.replace("\n","").split() for at in open(ini,"r").readlines() ]
	for at in range(len(fob)):
		try:
			if fob[at][0] == 'Index':
				if fob[at][1] == str(n) :
					geom = fob[at+2:at+2 + NAtoms]
					return geom
		except:
			pass

def write(dat,loc):
	fob1 = open(loc+"/geom","w+") 
	fob2 = open(loc+"/veloc","w+") 
	for i in dat:
		fob1.writelines(" ".join(i[:6]) +"\n" )
		if (zeroMomentum == False):
			fob2.writelines(" ".join(i[6:]) +"\n" ) 
		else:
			printline = np.zeros((3))
			fob2.writelines(" ".join(map(str, printline ) ) +"\n" )
	fob1.close()
	fob2.close()

dirs = []
NF = (nproc-off)
print (f"\tThere will be {NF} jobs.")

for i in range(off,nproc):
	name = f'TRAJ/traj-{i}'
	dirs.append(name)

#generate folders
dum = 0
deleteFlag = False
for i in range(off,nproc):
	if ( not os.path.exists(dirs[ dum ]) ):
		sp.call(f" mkdir -p {dirs[ dum ]} ", shell=True)
	else:
		if ( deleteFlag == False ):
			deleteQ = input("Found folders. Should I delete them all? 'Y' or 'N' ")
			if ( deleteQ == 'Y' ):
				print ("\tOkay. I will delete them.")
				deleteFlag = True
			else:
				print (f"Quitting since user gave {deleteQ} as answer.")
				exit()
	if ( deleteFlag == True ):	
		sp.call(f" rm -r {dirs[ dum ]} ", shell=True)
		sp.call(f" mkdir -p {dirs[ dum ]} ", shell=True)
	try:
		os.mkdir(dirs[ dum ] +"/QM" )
		os.mkdir(dirs[ dum ] +"/QM/temp" )
		os.mkdir(dirs[ dum ] +"/restart" )
	except :
		pass
	dum += 1

def makeinp(filename):
	I = open('input',"r").readlines()  
	for i in range(len(I)):
		if ( len(I[i].split()) > 0 and I[i].split()[0] == 'rngseed' ):
			I[i] = 'rngseed '+ str(random.randint(-1E8,1E8))   
		fob = open(filename,"w")
		fob.writelines(I) 
		fob.close()

cwd = os.getcwd()

#copy executables
dum = 0
for i in range(off,nproc):
	print (f"TRAJ: {i} of {nproc} ({dum+1})")
	makeinp(dirs[ dum ]+"/input") 
	#shutil.copy2('input',dirs[i-off])
	shutil.copy2('run.sh',dirs[ dum ])
	shutil.copy2('submit.sbatch',dirs[ dum ])
	#shutil.copy2('./QM/MOLPRO.resources',dirs[i]+"/QM/")
	resources = []
	for j in open("./QM/MOLPRO.resources","r").readlines():
		if (j.split()[0] == "scratchdir"):
			resources.append( "scratchdir " + cwd + "/%s/temp\n" %(dirs[ dum ]) )
		elif (j.split()[0] == "savedir"):
			resources.append( "savedir " + cwd + "/%s/restart\n" %(dirs[ dum ]) )
		elif (j.split()[0] == "memory"):
			resources.append( "memory %s\n" %memory )
		else:
			resources.append(j)
	open(dirs[ dum ]+"/QM/MOLPRO.resources","w+").writelines(resources)

	shutil.copy2('./QM/MOLPRO.template',dirs[ dum ]+"/QM/")
	shutil.copy2('./QM/runQM.sh',dirs[ dum ]+"/QM/") 
	#write(getxyz( dum ),dirs[ dum ]) # This one uses a different set of coords and velocs for each trajetory
	write(getxyz( i-off ),dirs[ dum ]) # This one uses the same coords and veloc for each set of traj-{} folders. This uses less initial conditions but is okay.
	os.chdir(dirs[ dum ])
	os.system("sbatch submit.sbatch")
	os.chdir(cwd)

	#	shutil.copy2('./QM/QM.in',dirs[i]+"/QM/")

	dum += 1
