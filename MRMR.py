#!/usr/bin/env python
import numpy as np
import os
import sys
import MDAnalysis as mda
import json

#subprocess.call(['python', 'script2.py', NameofSystem , deltatstring, tpr_Name, trajectory_Name, Project_Directory, proteintype ])
deltat=int(sys.argv[2])
tpr_Name=sys.argv[3]
trajectory_Name=sys.argv[4]
Project_Directory=sys.argv[5]
proteintype=sys.argv[6]
PDbfile=sys.argv[7]
Name=sys.argv[1].split()

rcut0=0.25
rcut1=0.6


#############################################################################################################################################
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
	vector=np.zeros((len(vector_0),len(vector_0)))
	for i in range(0,len(vector_0)):
		for j in range(0,len(vector_0)):
			vector[len(vector_0)-1-i][j]=vector_0[i][j]
	return vector

		
def UpperTriangleMatrix(vectorA, cutpoint0,cutpoint1):
	Matrix_2=[]  
	for i in range(0,len(vectorA)):
		Matrix_1=[]
		for j in range(0,len(vectorA)):
			if cutpoint0 <= float(vectorA[i][j]) <= cutpoint1:
				b=0
			else:
				b=1
			Matrix_1.append(b)
		Matrix_2.append(Matrix_1)	    	
	array_2=np.array(Matrix_2)
	inds = np.triu_indices_from(array_2, k = 1)
	data = array_2[inds].tolist()
	return data
		
def UpperTriangleIndex(vectorAsize,icomponent):
	vectorA=np.zeros((vectorAsize,vectorAsize))    	
	inds = np.triu_indices_from(vectorA, k = 1)
	a=inds[0][icomponent]
	b=inds[1][icomponent]
	return a,b
	
	
		
def vectorModify_protein(vector,Protein_type,Pr_Num):
	if len(vector) == int(Pr_Num):  
		rangev1=0
		rangev2=int(Pr_Num)		
	if len(vector) != int(Pr_Num):
		if Protein_type == "AA":
			rangev1=0
			rangev2=int(Pr_Num)
		if Protein_type == "BB":
			rangev1=int(Pr_Num)
			rangev2=2*int(Pr_Num)
		if Protein_type == "CC":
			rangev1=2*int(Pr_Num)
			rangev2=3*int(Pr_Num)
	print(rangev1,rangev2)
	vector_modify=np.zeros((int(Pr_Num),int(Pr_Num)))
	for i in range(rangev1,rangev2):
		for j in range(rangev1,rangev2):
			vector_modify[i-rangev1][j-rangev1]=vector[i][j]	
	return vector_modify

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

#########################################################################################################################################
def MIHI(className0,IndexNum0,IndexNum1):
	V1=[]
	V2=[]
	sum1=0
	#print(IndexNum0)
	#print(IndexNum1)	
	for j in range(0,len(className0[0].DistanceDist)):
		for ll in range(0,len(className0)):		
			V1.append(className0[ll].DistanceDist[j][IndexNum0])
			sum1=sum1+className0[ll].DistanceDist[j][IndexNum0]
			V2.append(className0[ll].DistanceDist[j][IndexNum1])
			sum1=sum1+className0[ll].DistanceDist[j][IndexNum1]	
	#print(V1)
	#print(V2)	
	if sum1 == 2*len(V1):
		kl_div_0=0		
	kl_div_0 = mi_fast(np.array(V1),np.array(V2)) 
	#print(kl_div_0)		
	return kl_div_0

	
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
	firstIDprotein=int(matrix[0][4])
	P_Num=int(matrix[len(matrix)-1][4])-firstIDprotein+1	
	vector2=[]
	for i in range(0,len(matrix)):
		vector2.append((matrix[i][4],matrix[i][3]))	
	for j in vector2:
		if int(j[0]) == IDNumer+firstIDprotein:
			bb=j[1] 
			break
	return bb

def PNumber(PDbfile):
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
	firstIDprotein=int(matrix[0][4])
	P_Num=int(matrix[len(matrix)-1][4])-firstIDprotein+1		
	return P_Num,firstIDprotein,matrix

#######################################################################################################################################################



P_Num=PNumber(PDbfile)[0]
print("P_Num",P_Num)
print("NAme",Name)

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
	

class MyClass:
    def __init__(self, class_attr, DistanceDist):
    	self.class_attr = class_attr
    	self.DistanceDist = DistanceDist

if proteintype == "Trimer" :
	Protein=["AA","BB","CC"]
else:
	Protein=["AA"]	
     
for pp in  Protein:
	className=['Class' + str(i) for i in range(len(Name))]        
	vectorName = ['vector' + str(i) for i in range(len(Name))] 
	vector_modify_Name = ['vector_modify' + str(i) for i in range(len(Name))]       
	conf_conf_dif_0_Name=['conf_conf_dif_' + str(i) for i in range(len(Name))]
	data_conf_conf_dif_Name=['data_conf_conf_dif_'+ str(i) for i in range(len(Name))]
	conf_conf_dif_0_Name=[[] for _ in range(len(Name))]
	
	
	for kk in range(steps):
		for ll in range(0,len(Name)):
			newdirectory2=os.path.join(Project_Directory,Name[ll])
			os.chdir(newdirectory2)
			vectorName[ll]=Data_fromFile("ContactMatrix/MatricContact"+str(kk)+".ndx")
			vector_modify_Name[ll]=vectorModify_protein(vectorName[ll],pp,P_Num)
			data_conf_conf_dif_Name[ll]=UpperTriangleMatrix(vector_modify_Name[ll],rcut0,rcut1)
			conf_conf_dif_0_Name[ll].append(data_conf_conf_dif_Name[ll])
								 	
	for classNumber in range(0,len(className)):
		className[classNumber]=MyClass(classNumber,conf_conf_dif_0_Name[classNumber])

	MI_values=[]		
	for i in range(0,len(className[0].DistanceDist[0])):
		V1=[]
		V2=[]
		for j in range(0,len(className[0].DistanceDist)):
			for classNum in range(0,len(Name)):
				V1.append(className[classNum].DistanceDist[j][i])
				V2.append(className[classNum].class_attr)
		MI_values.append((mi_fast(V1, V2),i))
	MI_sorted=sorted(MI_values, key=lambda x: (x[0]), reverse=True)[:100]
	print(MI_sorted)
	
	nn=int(MI_sorted[0][1])
	mm=float(MI_sorted[0][0])		
	MIMaxlist = MI_sorted
	MMR_value=[]
	while len(MIMaxlist)!=0:
		Value1=[]
		Value2=[]
		sumMIHI=0
		for Index_MMR in range(0,len(MIMaxlist)):
			if len(MMR_value) != 0:
				sumMIHI=MIHI(className,nn,int(MIMaxlist[Index_MMR][1]))
				MIHI_count=len(MMR_value)+1
				for item in range(0,len(MMR_value)):
					sumMIHI=sumMIHI+ MIHI(className,nn,int(MMR_value[item][1]))
										
			else:
				sumMIHI=MIHI(className,nn,int(MIMaxlist[Index_MMR][1]))
				MIHI_count=1				
			Value1.append(( +((MIHI(className,nn,int(MIMaxlist[Index_MMR][1])))-(sumMIHI/MIHI_count)),int(MIMaxlist[Index_MMR][1]) ))
			Value2=sorted(Value1, key=lambda x: (x[0]), reverse=True)
		MMR_value.append((Value2[0][0],Value2[0][1]))
		for sss in range(0,len(MIMaxlist)):
			if int(MIMaxlist[sss][1]) == Value2[0][1]:				
				MIMaxlist.pop(sss)
				break
		nn=Value2[0][1]
		mm=Value2[0][0]
	print("MMR_value:",MMR_value)	
	pair_Max=[]
	for pair in range(0,len(MMR_value)):		
		pair_Max.append((UpperTriangleIndex(len(vector_modify_Name[0]),MMR_value[pair][1])))
	print(pair_Max)
	print(len(pair_Max))
	
	os.chdir(Project_Directory) 
	VectorName=[]
	for j in range(0,len(pair_Max)):
		VectorName.append((ConvertIDtoName(pair_Max[j][0],PDbfile),ConvertIDtoName(pair_Max[j][1],PDbfile)))
	print(VectorName)
	print(len(VectorName))
	print("size",len(vectorName[0]))

	with open("MRMR_PAIRNAME_"+pp+".json", 'w', encoding='utf-8') as f:
		f.write('\n'.join('{} {}'.format(*tup) for tup in VectorName))
	
	with open("MRMR_PAIRNumber_"+pp+".json", 'w', encoding='utf-8') as f:
		f.write('\n'.join('{} {}'.format(*tup) for tup in pair_Max))
	
