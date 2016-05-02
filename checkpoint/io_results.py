import re
import os
import math
import numpy as np
import mdtraj as md
from msmbuilder.hmm import GaussianHMM as GHMM
from msmbuilder.decomposition import tICA
import subprocess as sub
from copy import deepcopy


########################################################################
#																	   #
#	 		   			   Print trajectories						   #
#																	   #
########################################################################


def print_trajs(trajectories, savepath):
#
	print("Printing trajectories...")
	
	wf = open(savepath, 'w')
	
	for i in range(len(trajectories)):
	#
		wf.write("Trajectory {:d} {:d}\n".format(len(trajectories[i]), len(trajectories[i][0])))
		
		for j in range(len(trajectories[i])):
		#
			wf.write("\tFrame\n".format(j))
			
			for k in range(len(trajectories[i][j])):
			#
				wf.write("\t\t{:f}\n".format(trajectories[i][j][k]))
			#
		#
	#
#


########################################################################
#																	   #
#	 		   			    Load trajectories						   #
#																	   #
########################################################################


def load_trajs(loadpath):
#
	trajectories = list()
	
	rf = open(loadpath, 'r')
	
	i = -1
	j = -1
	k = -1
	
	for line in rf:
	#
		m = re.match("\s\s([\-\d\.]+)", line)
		
		if m:
		#
			k += 1
					
			trajectories[i][j][k] = float(m.group(1))
			
			continue
		#
		else:
		#
			m = re.match("\sFrame", line)
			
			if m:
			#
				j += 1
				
				k = -1
				
				continue
			#
			else:
			#
				m = re.match("Trajectory\s(\d+)\s(\d+)", line)
				
				if m:
				#
					trajectories.append([[0.0 for j in range(int(m.group(2)))] for i in range(int(m.group(1)))])
					
					j = -1
					
					i += 1
					
					continue
				#
			#
		#
		
		raise Exception("Could not read line correctly.")
	#
	
	return trajectories
#


########################################################################
#																	   #
#	 		   			   Print studied features					   #
#																	   #
########################################################################


def print_feat(trajectory_atoms, savepath):
#
	print("Printing features...")
	
	wf = open(savepath, 'w')
	
	wf.write("\t".join(trajectory_atoms))
#


########################################################################
#																	   #
#	 		   			   Load studied features					   #
#																	   #
########################################################################


def load_feat(loadpath):
#
	rf = open(loadpath, 'r')
	
	features = None
	
	for line in rf:
	#
		features = re.split('\s+', line)
	#
	
	return features
#


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
	
	temp = deepcopy(tica_tot.eigenvectors_)
	
	evecs = np.asarray(temp).T.tolist() # .. # Writing eigenvectors as lines is more convenient in terms of data structure.
	
	for i in range(len(evecs)):
	#
		wf.write("\t\tevec {:d}\n".format(i))
		
		for j in range(len(evecs[i])):
		#
			wf.write("\t\t\t{:f}\n".format(evecs[i][j]))
		#
	#
	
	del temp
	del evecs
	
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
	
	tica_dict = dict()
	
	equil_values = {'True' : True, 'False' : False}
	
	attribute = None
	pos_i = -1
	pos_j = -1
	
	for line in rf:
	#
		m = re.match("\s\s\s([\-\d\.]+)", line)
		
		if m:
		#
			pos_j += 1
					
			tica_dict[attribute][pos_i][pos_j] = float(m.group(1))
			
			continue
		#
		else:
		#
			m = re.match("\s\s([\-\d\.]+)", line)
			
			if m:
			#
				pos_i += 1
				
				tica_dict[attribute][pos_i] = float(m.group(1))
				
				continue
			#
			else:
			#
				m = re.match("\s\s[\w\s]+", line)
				
				if m:
				#
					pos_i += 1
					
					pos_j = -1
					
					continue
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
			tica_dict['equil'] = equil_values[m.group(1)]
			
			continue
		#
		
		m = re.match("equil_dists : (\d+)", line)
		if m:
		#
			tica_dict['equil_dists'] = int(m.group(1))
			
			continue
		#
		
		m = re.match("usable_comps : (\d+)", line)
		if m:
		#
			tica_dict['usable_comps'] = int(m.group(1))
			
			continue
		#
		
		m = re.match("n_comp : (\d+)", line)
		if m:
		#
			tica_dict['n_comp'] = int(m.group(1))
			
			continue
		#
		
		m = re.match("frameskip : (\d+)", line)
		if m:
		#
			tica_dict['frameskip'] = int(m.group(1))
			
			continue
		#
		
		m = re.match("n_features : (\d+)", line)
		if m:
		#
			tica_dict['n_feat'] = int(m.group(1))
			
			continue
		#
		
		m = re.match("tICA", line)
		if m:
		#
			tica_dict['eigenvalues'] = [0.0 for i in range(tica_dict['n_feat'])]
			tica_dict['eigenvectors'] = [[0.0 for i in range(tica_dict['n_feat'])] for i in range(tica_dict['n_feat'])] # .. # REMEMBER : This is the transpose; its eigenvectors are lines.
			tica_dict['means'] = [0.0 for i in range(tica_dict['n_feat'])]
			tica_dict['timescales'] = [0.0 for i in range(tica_dict['n_feat'])]
			tica_dict['n_observations'] = [0.0]
			tica_dict['n_sequences'] = [0.0]
				
			continue
		#
		
		m = re.match("$", line)
		if m:
		#
			continue
		#
		
		raise Exception("Could not read line correctly.")
	#
	
	return tica_dict
#


########################################################################
#																	   #
#	 		   			  Print tICA projections					   #
#																	   #
########################################################################


def print_tica_projs(ic_projs, n_comp, savepath):
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
			
			for k in range(n_comp):
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
#	 		   		 Print the results of a tICA scan				   #
#																	   #
########################################################################


def print_scan(compares, centers, out_root):
#
	for i in range(len(compares)):
	#
		with open("{:s}_scan_EVec_{:d}.txt".format(out_root, i+1), 'w') as wf:
		#
			wf.write("X:Centers\n")
			
			for center in centers:
			#
				wf.write("{:d}\n".format(center))
			#
			
			wf.write("Cos_angles\n")
			
			for cos in compares[i][0]:
			#
				wf.write("{:f}\n".format(math.fabs(cos)))
			#
			
			wf.write("R_t\n")
			
			for rt in compares[i][1]:
			#
				wf.write("{:f}\n".format(rt))
			#
		#
	#
#


########################################################################
#																	   #
#	 		   			  	Print HMM results						   #
#																	   #
########################################################################


def print_hmm(hmm, savepath):
#
	wf = open(savepath, 'w')
	
	wf.write("n_states : {:d}\n".format(len(hmm.means_)))
	wf.write("n_dims : {:d}\n".format(len(hmm.means_[0])))
	
	wf.write("Log-likelihood : {:f}\n".format(hmm.fit_logprob_[-1]))
	
	wf.write("Means :\n")
	
	for i in range(len(hmm.means_)):
	#
		wf.write("\tMean {:d}\n".format(i))
		
		for j in range(len(hmm.means_[i])):
		#
			wf.write("\t\t{:f}\n".format(hmm.means_[i][j]))
		#
	#
	
	wf.write("Vars :\n")
	
	for i in range(len(hmm.vars_)):
	#
		wf.write("\tVar {:d}\n".format(i))
		
		for j in range(len(hmm.vars_[i])):
		#
			wf.write("\t\t{:f}\n".format(hmm.vars_[i][j]))
		#
	#
	
	wf.write("Populations :\n")
	
	for i in range(len(hmm.populations_)):
	#
		wf.write("\t{:f}\n".format(hmm.populations_[i]))
	#
	
	wf.write("Transmat :\n")
	
	for i in range(len(hmm.transmat_)):
	#
		wf.write("\tState {:d}\n".format(i))
		
		for j in range(len(hmm.transmat_[i])):
		#
			wf.write("\t\t{:f}\n".format(hmm.transmat_[i][j]))
		#
	#
	
	wf.write("Timescales :\n")
	
	for i in range(len(hmm.timescales_)):
	#
		wf.write("\t{:f}\n".format(hmm.timescales_[i]))
	#
#


########################################################################
#																	   #
#	 		   			  	Load HMM results						   #
#																	   #
########################################################################


def load_hmm(loadpath):
#
	rf = open(loadpath, 'r')
	
	hmm_dict = dict()
	
	attribute = str()
	pos_i = -1
	pos_j = -1
	
	for line in rf:
	#
		print(line)
		
		m = re.match("\s\s([-\d\.]+)", line)
		
		if m:
		#
			pos_j += 1
					
			hmm_dict[attribute][pos_i][pos_j] = float(m.group(1))
			
			continue
		#
		else:
		#
			m = re.match("\s([-\d\.]+)", line)
			
			if m:
			#
				pos_i += 1
				
				hmm_dict[attribute][pos_i] = float(m.group(1))
				
				continue
			#
			else:
			#
				m = re.match("\s[\w\s]+", line)
				
				if m:
				#
					pos_i += 1
					
					pos_j = -1
					
					continue
				#
				else:
				#
					pos_i = -1
					
					m = re.match("n_states : (\d+)", line)
					if m:
					#
						#print("n_states")
						
						hmm_dict['n_states'] = int(m.group(1))
						
						continue
					#
					
					m = re.match("n_dims : (\d+)", line)
					if m:
					#
						#print("n_dims")
						
						hmm_dict['n_dims'] = int(m.group(1))
						
						continue
					#
					
					m = re.match("Log-likelihood : ([-\d\.]+)", line)
					if m:
					#
						#print("logL")
						
						hmm_dict['logL'] = float(m.group(1))
						
						continue
					#
					
					m = re.match("Means :", line)
					if m:
					#
						#print("means")
						
						hmm_dict['means'] = [[0.0 for i in range(hmm_dict['n_dims'])] for j in range(hmm_dict['n_states'])]
						
						attribute = 'means'
						
						continue
					#
					
					m = re.match("Vars :", line)
					if m:
					#
						#print("vars")
						
						hmm_dict['vars'] = [[0.0 for i in range(hmm_dict['n_dims'])] for j in range(hmm_dict['n_states'])]
						
						attribute = 'vars'
						
						continue
					#
					
					m = re.match("Populations :", line)
					if m:
					#
						#print("popls")
						
						hmm_dict['popls'] = [0.0 for i in range(hmm_dict['n_states'])]
						
						attribute = 'popls'
						
						continue
					#
					
					m = re.match("Timescales :", line)
					if m:
					#
						#print("timescales")
						
						hmm_dict['timescales'] = [0.0 for i in range(hmm_dict['n_states'])]
						
						attribute = 'popls'
						
						continue
					#
					
					m = re.match("Transmat :", line)
					if m:
					#
						#print("transmat")
						
						hmm_dict['transmat'] = [[0.0 for i in range(hmm_dict['n_states'])] for j in range(hmm_dict['n_states'])]
						
						attribute = 'transmat'
						
						continue
					#
				#
			#
		#
		
		raise Exception("Could not read line correctly.")
	#
	
	return hmm_dict
#


########################################################################
#																	   #
#	 		   		     	 Print centroids						   #
#																	   #
########################################################################


def print_centroids(hmm, trajins, structin, timeskip, state_centers, trajout):
#
	print("Printing centroids...")
	
	with open("{:s}_states.dat".format(trajout), 'w') as wf:
	#
		for i in range(len(state_centers)):
		#
			if structin is None:
			#
				with open(trajins[state_centers[i][0][0]], 'r') as rf:
				#
					with open("{:s}_state_{:d}.pdb".format(trajout, i+1), 'w') as sf:
					#
						trigger = False
			
						for line in rf:
						#
							if trigger:
							#
								sf.write(line)
					
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
									sf.write(line)
									
									wf.write("\n\nState {:d}\n".format(i+1))
									
									wf.write("State mean : \n")
									wf.write(", ".join(map(str, hmm.means_[i])))
									wf.write("\n")
						
									wf.write("State variance : \n")
									wf.write(", ".join(map(str, hmm.vars_[i])))
									wf.write("\n")
						
									wf.write("State population : \n")
									wf.write(str(hmm.populations_[i]))
									wf.write("\n")
						
									trigger = True
								#
							#
						#
					#
				#
			#
			else:
			#
				wf.write("\n\nGROMACS command for state {:d}\n".format(i+1))
				wf.write("gmx trjconv -f {:s} -s {:s} -o {:s}_state_{:d}.pdb -pbc mol -ur compact -b {:f} -e {:f}\n\n".format(trajins[state_centers[i][0][0]], structin, trajout, i+1, (float(state_centers[i][0][1]) - 0.5)*timeskip, (float(state_centers[i][0][1]) + 0.5)*timeskip))
			
				#wf.open("{:s}_states.dat".format(trajout), 'a')
			
				#wf.write("\nState {:d}\n\n".format(i+1))
				wf.write("State {:d}\n".format(i+1))
			
				#wf.write("State mean : ")
				#wf.write(", ".join(map(str, hmm.means_[i])))
				#wf.write("\n")
				wf.write("State mean : \n")
				wf.write(", ".join(map(str, hmm.means_[i])))
				wf.write("\n")
			
				#wf.write("State variance : ")
				#wf.write(", ".join(map(str, hmm.vars_[i])))
				#wf.write("\n")
				wf.write("State variance : \n")
				wf.write(", ".join(map(str, hmm.vars_[i])))
				wf.write("\n")
			
				#wf.write("State population : ")
				#wf.write(str(hmm.populations_[i]))
				#wf.write("\n")
				wf.write("State population : \n")
				wf.write(str(hmm.populations_[i]))
				wf.write("\n")
			
				#wf.close()
			#
		#
	#
#
