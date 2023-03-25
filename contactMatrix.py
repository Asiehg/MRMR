#!/usr/bin/env python
#!/usr/bin/env python
import os
import sys
import MDAnalysis as mda

deltat=1000 #ps
tpr_Name="step7_production.tpr"
trajectory_Name="step7_production.trr"

# Load trajectory and topology files
u = mda.Universe(tpr_Name, trajectory_Name)

# Get the number of frames and time step
num_frames = len(u.trajectory)
time_step = u.trajectory.dt

# Calculate the last frame time
last_frame_time = (num_frames - 1) * time_step
steps=int((last_frame_time/deltat) +1)


def unquote(s):
    return s[1+s.find('"'):s.rfind('"')]

def uncomment(s):
    return s[2+s.find('/*'):s.rfind('*/')]

def col(c):
    color = c.split('/*')
    value = unquote(color[1])
    color = unquote(color[0]).split()
    sys.stderr.write("%s: %s %s %s\n"%(c.strip(),color[0], color[1], value))
    return color[0], value

#print(sys.argv[1])
# Open the xpm file for reading

def Convertxpm(ivalue):
	xpm = open("distframes"+str(ivalue)+".xpm")
	meta = [xpm.readline()]
	while not meta[-1].startswith("static char *gromacs_xpm[]"):
	    meta.append(xpm.readline())

	dim = xpm.readline()
	nx, ny, nc, nb = [int(i) for i in unquote(dim).split()]
	colors = dict([col(xpm.readline()) for i in range(nc)])
	f = open("ContactMatrix/MatricContact"+str(ivalue)+".ndx", "w")
	for i in xpm:
	    if i.startswith("/*"):
	    	continue 
	    j = unquote(i)    
	    z = [colors[j[k:k+nb]] for k in range(0,nx*nb,nb)]
	    f.write(" ".join(z)+"\n")
	    
for i in range(steps):
	timeframe = deltat * i
	os.system('gmx_mpi select -s %s -on protein -select "group Protein"' %(tpr_Name))
	string1="distframes"+str(i)+""
	os.system('gmx_mpi mdmat -s %s -f %s -n protein -frames %s -t 7 -b %s -e %s' %(tpr_Name,trajectory_Name,string1,timeframe,timeframe))
	Convertxpm(i)

