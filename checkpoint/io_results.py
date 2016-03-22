import re
import sys
import math
import numpy as np
import mdtraj as md
from msmbuilder.hmm import GaussianHMM as GHMM
from msmbuilder.decomposition import tICA


########################################################################
#																	   #
#	 		   				Print tICA results						   #
#																	   #
########################################################################


def print_tica(tica_tot, equil, equil_dists, n_comp, usable_comps, frameskip, savepath):
#
	print("Printing tICA...")
	
	wf = open(savepath, 'w')
	
	wf.write("equil : ")
	if equil:
	#
		wf.write("True\n")
	#
	else:
	#
		wf.write("False\n")
	#
	
	wf.write("equil_dists : {:d}\n".format(equil_dists))
	wf.write("usable_comps : {:d}\n".format(usable_comps))
	wf.write("n_comp : {:d}\n".format(n_comp))
	wf.write("frameskip : {:d}\n".format(frameskip))
	wf.write("n_features : {:d}\n".format(len(tica_tot.means_)))
	wf.write("tICA\n")
	
	wf.write("\teigenvalues\n")
	for i in range(len(tica_tot.eigenvalues_)):
	#
		wf.write("\t\t{:f}\n".format(tica_tot.eigenvalues_[i]))
	#
	
	wf.write("\teigenvectors\n")
	for i in range(len(tica_tot.eigenvectors_)):
	#
		wf.write("\t\tevec {:d}\n".format(i))
		
		for j in range(len(tica_tot.eigenvectors_[i])):
		#
			wf.write("\t\t\t{:f}\n".format(tica_tot.eigenvectors_[i][j]))
		#
	#
	
	wf.write("\tmeans\n")
	for i in range(len(tica_tot.means_)):
	#
		wf.write("\t\t{:f}\n".format(tica_tot.means_[i]))
	#
	
	wf.write("\tn_observations\n\t\t{:d}\n".format(tica_tot.n_observations_))
	wf.write("\tn_sequences\n\t\t{:d}\n".format(tica_tot.n_sequences_))
	
	wf.write("\ttimescales\n")
	for i in range(len(tica_tot.timescales_)):
	#
		if tica_tot.timescales_[i] != tica_tot.timescales_[i]:
		#
			wf.write("\t\t0.0\n")
		#
		else:
		#
			wf.write("\t\t{:f}\n".format(tica_tot.timescales_[i]))
		#
	#
#


########################################################################
#																	   #
#	 		   				Load tICA results						   #
#																	   #
########################################################################


def load_tica(loadpath):
#
	rf = open(loadpath, 'r')
	
	side_dict = {'equil' : None, 'equil_dists' : 0, 'n_comp' : 0, 'usable_comps' : 0, 'frameskip' : 0, 'n_features' : 0}
	
	tica_tot = None
	tica_dict = None
	
	equil_values = {'True' : True, 'False' : False}
	
	attribute = None
	pos_i = -1
	pos_j = -1
	
	for line in rf:
	#
		m = re.match("\s\s\s([\d\.]+)", line)
		
		if m:
		#
			pos_j += 1
					
			tica_dict[attribute][pos_i][pos_j] = float(m.group(1))
		#
		else:
		#
			m = re.match("\s\s([\d\.]+)", line)
			
			if m:
			#
				pos_i += 1
				
				tica_dict[attribute][pos_i] = float(m.group(1))
			#
			else:
			#
				m = re.match("\s\s[\w\s]+", line)
				
				if m:
				#
					pos_i += 1
					
					pos_j = -1
				#
				else:
				#
					m = re.match("\s([\w\_]+)", line)
					
					if m:
					#
						attribute = m.group(1)
						
						pos_i = -1
						
						continue
					#
				#
			#
		#
		
		m = re.match("equil : (\w+)", line)
		if m:
		#
			side_dict['equil'] = equil_values[m.group(1)]
			
			continue
		#
		
		m = re.match("equil_dists : (\d+)", line)
		if m:
		#
			side_dict['equil_dists'] = int(m.group(1))
			
			continue
		#
		
		m = re.match("usable_comps : (\d+)", line)
		if m:
		#
			side_dict['usable_comps'] = int(m.group(1))
			
			continue
		#
		
		m = re.match("n_comp : (\d+)", line)
		if m:
		#
			side_dict['n_comp'] = int(m.group(1))
			
			continue
		#
		
		m = re.match("frameskip : (\d+)", line)
		if m:
		#
			frameskip = int(m.group(1))
			
			continue
		#
		
		m = re.match("n_features : (\d+)", line)
		if m:
		#
			n_features = int(m.group(1))
			
			continue
		#
		
		m = re.match("tICA", line)
		if m:
		#
			if side_dict['equil'] is None or side_dict['equil_dists'] == 0 or side_dict['usable_comps'] == 0 or side_dict['n_comp'] == 0 or side_dict['frameskip'] == 0 or side_dict['n_feat'] == 0:
			#
				raise Exception("tICA has not been printed/read correctly")
			#
			else:
			#
				tica_dict = {'eigenvalues' : [0.0 for i in range(side_dict['n_feat'])], 'eigenvectors' : [[0.0 for i in range(side_dict['n_feat'])] for i in range(side_dict['n_feat'])], 'means' : [0.0 for i in range(side_dict['n_feat'])], 'timescales' : [0.0 for i in range(side_dict['n_feat'])], 'n_observations' : [0.0], 'n_sequences' : [0.0]}
				
				continue
			#
		#
		
		raise Exception("Could not read line correctly.")
	#
	
	return side_dict, tica_dict
#


########################################################################
#																	   #
#	 		   			  Print tICA projections					   #
#																	   #
########################################################################


def print_tica_projs(ic_projs, savepath):
#
	print("Printing projections on tICA components...")
	
	wf = open(savepath, 'w')
	
	wf.write("Number of trajectories : {:d}\n".format(len(ic_projs)))
	
	for i in range(len(ic_projs)):
	#
		wf.write("Trajectory {:d}\n".format(i+1))
		
		for j in range(len(ic_projs[i])):
		#
			wf.write("\tFrame {:d}\n".format(j+1))
			
			for k in range(len(ic_projs[i][j])):
			#
				wf.write("\t\t{:f}\n".format(ic_projs[i][j][k]))
			#
		#
	#
#

#def load_tica_projs(loadpath):
#
	
	
	#return ic_splits
#


########################################################################
#																	   #
#	 		   			  	Print HMM results						   #
#																	   #
########################################################################





########################################################################
#																	   #
#	 		   		     	 Print centroids						   #
#																	   #
########################################################################


def print_centroids(hmm, trajins, topin, timeskip, state_centers, trajout):
#
	print("Printing centroids...")
	
	for i in range(len(state_centers)):
	#
		if topin is None:
		#
			rf = open(trajins[state_centers[i][0][0]], 'r')
			
			wf = open("{:s}_state_{:d}.pdb".format(trajout, i+1), 'w')
			
			trigger = False
			
			for line in rf:
			#
				if trigger:
				#
					wf.write(line)
					
					if re.match("ENDMDL", line):
					#
						trigger = False
					#
				#
				else:
				#
					m = re.match("MODEL\s+(\d+)", line)
					
					if m and int(m.group(1)) == state_centers[i][0][1] + 1:
					#
						wf.write(line)
						
						wf.write("REMARK   6 State mean : ")
						wf.write(", ".join(map(str, hmm.means_[i])))
						wf.write("\n")
						
						wf.write("REMARK   6 State variance : ")
						wf.write(", ".join(map(str, hmm.vars_[i])))
						wf.write("\n")
						
						wf.write("REMARK   6 State population : ")
						wf.write(str(hmm.populations_[i]))
						wf.write("\n")
						
						trigger = True
					#
				#
			#
		#
	#
#
