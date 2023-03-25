#!/usr/bin/env python
import numpy as np
import os
import sys
import MDAnalysis as mda
import json

deltat=int(sys.argv[2])
tpr_Name=sys.argv[3]
trajectory_Name=sys.argv[4]
Project_Directory=sys.argv[5]
PDbfile=sys.argv[6]
Name=sys.argv[1].split()

#------------------------------------------------------------------------------------------------------------
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




def mi_fast(x, y):
    n_samples = len(x)
    bin_range = [min(np.min(x), np.min(y)), max(np.max(x), np.max(y))]
    n_bins = bin_range[1] - bin_range[0]+1
    xy_c = np.histogram2d(x, y, bins=n_bins, range=[bin_range, bin_range])[0]
    x_c = np.sum(xy_c, 1)
    y_c = np.sum(xy_c, 0)
    Hxy = entropy(np.ravel(xy_c)) / n_samples + np.log(n_samples)
    Hx = entropy(np.ravel(x_c)) / n_samples + np.log(n_samples)
    Hy = entropy(np.ravel(y_c)) / n_samples + np.log(n_samples)    
    return Hx + Hy - Hxy

def entropy(arr):
    return np.sum([-p * np.log(p) if p else 0 for p in arr]) 
    
def Data_fromFile(FolderName):
	f1 = open(FolderName,"r")
	lines1 = f1.readlines()
	linesarr1= np.array(lines1)
	f1.close()
	vector_0=[]
	for i in range(0,len(lines1)):
		linevalue=linesarr1[i].split(" ")
		resvalue = [ele for ele in linevalue if ele != '']
		vector_0.append(resvalue)
	size_vector0=len(vector_0) 
	vector=np.zeros((len(vector_0),len(vector_0)))
	for ii in range(0,len(vector_0)):
		for jj in range(0,len(vector_0)):
			vector[len(vector_0)-1-ii][jj]=vector_0[ii][jj]			
	inds = np.triu_indices_from(vector, k = 1)
	data = vector[inds].tolist()
	return data,size_vector0

def BiningDistance(vectorA):
	VectorB=[]
	for i in range(0,len(vectorA)):
		M=np.arange(0,7.1,0.1)
		for mem in range(0,len(M)-1):
			if M[mem] <= float(vectorA[i]) <= M[mem+1]:
				VectorB.append(int(10*M[mem+1]))													
	return VectorB	

def UpperTriangleIndex(vectorAsize,icomponent):
	vectorA=np.zeros((vectorAsize,vectorAsize))    	
	inds = np.triu_indices_from(vectorA, k = 1)
	a=inds[0][icomponent]
	b=inds[1][icomponent]
	return a,b
	
def ConvertIDtoName(IDNumer,PDbfile):
	file = open(PDbfile, "r")
	prev_line = ""
	data = []
	while True:
	    line = file.readline()
	    if not line:
	    	break
	    if line.startswith("ATOM"):
	    	row = [x for x in line.split(" ") if x != '']
	    	data.append(row)

	    if line.startswith("TER"):
	    	break
	    row = [x for x in line.split(" ") if x != '']
	    data.append(row)
	matrix = np.array(data)
	vector2=[]
	for i in range(0,len(matrix)):
		vector2.append((matrix[i][4],matrix[i][3]))
	firstIDprotein=int(matrix[0][4])		
	for j in vector2:
		if int(j[0]) == IDNumer+firstIDprotein:
			bb=j[1] 
			break
	return bb
	


def ContactMatrix(Name,steps,tpr_Name,trajectory_Name,deltat):
	# specify the directory name you want to create
	directory = "ContactMatrix"
	# use try-except block to handle errors
	try:
	    # create the directory
	    os.mkdir(directory)
	    print(f"Directory '{directory}' created successfully!")
	except OSError as error:
	    print(f"Directory '{directory}' already exists or could not be created: {error}")
	os.system('gmx_mpi select -s %s -on protein -select "group Protein"' %(tpr_Name))	
	for i in range(steps):
		timeframe = deltat * i		
		string1="distframes"+str(i)+""
		os.system('gmx_mpi mdmat -s %s -f %s -n protein -frames %s -t 7 -b %s -e %s' %(tpr_Name,trajectory_Name,string1,timeframe,timeframe))
		Convertxpm(i)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------	

for NameItem in Name:
	newdirectory=os.path.join(Project_Directory,NameItem)
	os.chdir(newdirectory)
	# Load trajectory and topology files
	u = mda.Universe(tpr_Name, trajectory_Name)

	# Get the number of frames and time step
	num_frames = len(u.trajectory)
	time_step = u.trajectory.dt

	# Calculate the last frame time
	last_frame_time = (num_frames - 1) * time_step
	steps=int((last_frame_time/deltat) +1)
	ContactMatrix(NameItem,steps,tpr_Name,trajectory_Name,deltat)


  	
vector= ['vector'+str(i) for i in range(len(Name))]
size_vector=['size_vector'+str(i) for i in range(len(Name))] 
vectorBinned= ['vectorBinned'+str(i) for i in range(len(Name))]  
 	
Data1=[]
for kk in range(steps):
	for ll in range(0,len(Name)): 
		newdirectory=os.path.join(Project_Directory,Name[ll])
		os.chdir(newdirectory)  	 	
		vector[ll],size_vector[ll]=Data_fromFile("ContactMatrix/MatricContact"+str(kk)+".ndx")
		vectorBinned[ll]=BiningDistance(vector[ll])
		Data1.append((vectorBinned[ll],ll))

os.chdir(Project_Directory) 
MI_values=[]
for jj in range(0,len(Data1[0][0])):
	v1=[]
	v2=[]
	for ii in range(0,len(Data1)):
		v1.append(Data1[ii][0][jj])
		v2.append(Data1[ii][1])
	MI_values.append((mi_fast(v1, v2),jj))

MI_sorted=sorted(MI_values, key=lambda x: (x[0]), reverse=True)[:100]

pair_Max=[]
for pair in range(0,len(MI_sorted)):		
	pair_Max.append((UpperTriangleIndex(size_vector[0],MI_sorted[pair][1])))		
VectorName=[]
for j in range(0,len(pair_Max)):
	VectorName.append((ConvertIDtoName(pair_Max[j][0],PDbfile),ConvertIDtoName(pair_Max[j][1],PDbfile)))	


RangeMi=int(len(Name)*(len(Name)-1)/2)
MI_valuess= ['MI_values_'+str(i) for i in range(0,RangeMi)] 
MI_valuess=[[] for _ in range(0,RangeMi)] 
MI_sorted_IND= ['MI_sorted_IND_'+str(i) for i in range(0,RangeMi)]

MIItem=-1
if len(Name)>2:
	for ll in range(0,len(Name)): 
		for ss in range(0,len(Name)): 
			if ll < ss: # skip same comparison
				MIItem=int(MIItem+1)
				for jj in range(0,len(Data1[0][0])):
					v1=[]
					v2=[]
					for ii in range(0,len(Data1)):
						if Data1[ii][1] == ll or Data1[ii][1] == ss:
							v1.append(Data1[ii][0][jj])
							v2.append(Data1[ii][1])
					MI_valuess[MIItem].append((mi_fast(v1, v2),jj))
						
							
for MIItem in range(0,len(MI_valuess)):
	MI_sorted_IND[MIItem]=sorted(MI_valuess[MIItem], key=lambda x: (x[0]), reverse=True)[:100]

with open("MI_All.json", 'w', encoding='utf-8') as f:
	f.write('\n'.join('{} {}'.format(*tup) for tup in MI_sorted))

with open("MI_PAIRNAME_All.json", 'w', encoding='utf-8') as f:
	f.write('\n'.join('{} {}'.format(*tup) for tup in VectorName))
	
with open("MI_PAIRNumber_All.json", 'w', encoding='utf-8') as f:
	f.write('\n'.join('{} {}'.format(*tup) for tup in pair_Max))

pair_Maxs = ['pair_Maxs_' + str(i) for i in range(0,RangeMi)]
pair_Maxs=[[] for _ in range(0,RangeMi)] 
VectorNames = ['VectorName_' + str(i) for i in range(0,RangeMi)]
VectorNames=[[] for _ in range(0,RangeMi)] 


for MIItem in range(0,RangeMi):
	for pair in range(0,len(MI_sorted_IND[0])):		
		pair_Maxs[MIItem].append((UpperTriangleIndex(size_vector[0],MI_sorted_IND[MIItem][pair][1])))
		
for MIItem in range(0,RangeMi):
	for j in range(0,len(pair_Maxs[MIItem])):
		VectorNames[MIItem].append((ConvertIDtoName((pair_Maxs[MIItem][j][0]),PDbfile),ConvertIDtoName((pair_Maxs[MIItem][j][1]),PDbfile)))	
	
MIItem=-1
for ll in range(0,len(Name)): 
	for ss in range(0,len(Name)):
		if ll < ss: # skip same comparison 
			MIItem=int(MIItem+1)
			with open("MI_Between_"+Name[ll]+"_and_"+Name[ss]+".json", 'w', encoding='utf-8') as f:
				f.write('\n'.join('{} {}'.format(*tup) for tup in MI_sorted_IND[MIItem]))
			with open("MI_PAIRNAME_Between_"+Name[ll]+"_and_"+Name[ss]+"-v5.json", 'w', encoding='utf-8') as f:
				f.write('\n'.join('{} {}'.format(*tup) for tup in VectorNames[MIItem]))
			with open("MI_PAIRNumber_Between_"+Name[ll]+"_and_"+Name[ss]+"-v5.json", 'w', encoding='utf-8') as f:
				f.write('\n'.join('{} {}'.format(*tup) for tup in pair_Maxs[MIItem]))			


