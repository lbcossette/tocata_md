import re
import sys
import math
import numpy as np
import mdtraj as md
from msmbuilder.hmm import GaussianHMM as GHMM


########################################################################
#																	   #
#	 		   				Print tICA results						   #
#																	   #
########################################################################


def print_tica(tica_tot, equil, equil_dists, n_comp, usable_comps, frameskip, savepath):
#
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
	
	wf.write("\tcomponents\n")
	
	for i in range(len(tica_tot.components_)):
	#
		wf.write("\t\tcomp {:d}\n".format(i))
		
		for j in range(len(tica_tot.components_[i]))
		#
			wf.write("\t\t\t{:f}\n".format(tica_tot.components_[i][j]))
		#
	#
	
	wf.write("\toffset_correlation\n")
	
	for i in range(len(tica_tot.offset_correlation_)):
	#
		wf.write("\t\trow {:d}\n".format(i))
		
		for j in range(len(tica_tot.offset_correlation_[i]))
		#
			wf.write("\t\t\t{:f}\n".format(tica_tot.offset_correlation_[i][j]))
		#
	#
	
	wf.write("\teigenvalues\n")
	
	for i in range(len(tica_tot.eigenvalues_)):
	#
		wf.write("\t\t{:f}\n".format(tica_tot.eigenvalues_[i]))
	#
	
	wf.write("\teigenvectors\n")
	
	for i in range(len(tica_tot.eigenvectors_)):
	#
		wf.write("\t\tevec {:d}\n".format(i))
		
		for j in range(len(tica_tot.eigenvectors_[i]))
		#
			wf.write("\t\t\t{:f}\n".format(tica_tot.eigenvectors_[i][j]))
		#
	#
	
	wf.write("\tmeans\n")
	
	for i in range(len(tica_tot.means_)):
	#
		wf.write("\t\t{:f}\n".format(tica_tot.means_[i]))
	#
	
	wf.write("\tn_observations : {:d}\n".format(tica_tot.n_observations_))
	
	wf.write("\tn_sequences : {:d}\n".format(tica_tot.n_sequences_))
	
	wf.write("\ttimescales\n")
	
	for i in range(len(tica_tot.timescales_)):
	#
		wf.write("\t\t{:f}\n".format(tica_tot.timescales_[i]))
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
		
		m = re.match("equil_dists : (\d+)", line)
		
		if m:
		#
			equil_dists = int(m.group(1))
			
			continue
		#
		
		m = re.match("usable_comps : (\d+)", line)
		
		if m:
		#
			usable_comps = int(m.group(1))
			
			continue
		#
		
		m = re.match("n_comp : (\d+)", line)
		
		if m:
		#
			n_comp = int(m.group(1))
			
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
			if equil is None or equil_dists is None or usable_comps is None or n_comp is None or frameskip is None or n_feat is None
			#
				raise Exception("tICA has not been printed correctly")
			#
			else:
			#
				tica_dict = {'components' : [[0.0 for i in range(n_feat)] for i in range(n_feat)], 'offset_correlation' : [[0.0 for i in range(n_feat)] for i in range(n_feat)], 'eigenvalues' : [0.0 for i in range(n_feat)], 'eigenvectors' : [[0.0 for i in range(n_feat)] for i in range(n_feat)], 'means' : [0.0 for i in range(n_feat)], 'timescales' : [0.0 for i in range(n_feat)], 'n_observations' : 0, 'n_sequences' : 0}
				
				continue
			#
		#
		
		raise Exception("Could not read line correctly.")
	#
	
	return tica_tot, equil, equil_dists, n_comp, usable_comps, frameskip
#


########################################################################
#																	   #
#	 		   			  Print tICA projections					   #
#																	   #
########################################################################


def print_tica_projs(ic_splits, savepath)
#
	
#

def load_tica_projs(loadpath)
#
	
	
	return ic_splits
#


########################################################################
#																	   #
#	 		   			  	Print HMM results						   #
#																	   #
########################################################################





########################################################################
#																	   #
#	 		   		     Print HMM fusion results					   #
#																	   #
########################################################################



