# script1.py
import subprocess


deltat=20000                                  #It is calculated contact matrix for every 20000ps
tpr_Name="step7_production.tpr"               #name of tpr of the systems
trajectory_Name="step7_production.trr"        #name of trr of the systems
Project_Directory="/scratch/fps7nd/MIFILES/"  #main directory that project exist                     
NameofSystem="BRER BtuB_POPC LPS_BtuB"        #Name of the systems that you wish to calculate MI
pdbfile ="step5_input.pdb"                    #Name of pdb file


deltatstring=str(deltat)
subprocess.call(['python', 'MI.py', NameofSystem , deltatstring, tpr_Name, trajectory_Name, Project_Directory , pdbfile])
