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
	
	wf.write("In equilibrium : ")
	
	if equil:
	#
		wf.write("True\n")
	#
	else:
	#
		wf.write("False\n")
	#
	
	wf.write("Equilibrium distributions : {:d}\n".format(equil_dists))
	
	wf.write("Amount of usable components : {:d}\n".format(usable_comps))
	
	wf.write("Amount of relevant components : {:d}\n".format(n_comp))
	
	wf.write("Number of frames skipped : {:d}\n".format(frameskip))
	
	wf.write("Number of features : {:d}\n".format(len(tica_tot.means_)))
	
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
	
	equil = None 
	equil_dists = None
	n_comp = None
	usable_comps = None 
	frameskip = None
	n_feat = None
	tica_tot = None
	
	components = None
	offset_correlation = None
	eigenvalues = None
	eigenvectors = None
	means = None
	n_observations = None
	n_sequences = None
	timescales = None
	
	attribute = None
	pos_i = 0
	pos_j = 0
	
	for line in rf:
	#
		m = re.match("\s", line)
		
		if m:
		#
			m = re.match("\s\s", line)
			
			if m:
			#
				m = re.match("\s\s\s", line)
				
				if m:
				#
					
				#
				else:
				#
					
				#
			#
			else:
			#
				m = re.match("\s([\w\_]+)", line)
				
				if m:
				#
					attribute = m.group(1)
				#
			#
		#
		
		m = re.match("In equilibrium : ", line)
		
		if m:
		#
			if re.search("True", line):
			#
				equil = True
			#
			elif re.search("False", line):
			#
				equil = False
			#
			else:
			#
				raise Exception("Equilibrium value could not be parsed.")
			#
			
			continue
		#
		
		m = re.match("Equilibrium distributions : (\d+)", line)
		
		if m:
		#
			equil_dists = int(m.group(1))
			
			continue
		#
		
		m = re.match("Amount of usable components : (\d+)", line)
		
		if m:
		#
			usable_comps = int(m.group(1))
			
			continue
		#
		
		m = re.match("Amount of relevant components : (\d+)", line)
		
		if m:
		#
			n_comp = int(m.group(1))
			
			continue
		#
		
		m = re.match("Number of frames skipped : (\d+)", line)
		
		if m:
		#
			frameskip = int(m.group(1))
			
			continue
		#
		
		m = re.match("Number of features : (\d+)", line)
		
		if m:
		#
			n_feat = int(m.group(1))
			
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
				components = [[0.0 for i in range(n_feat)] for i in range(n_feat)]
				offset_correlation = [[0.0 for i in range(n_feat)] for i in range(n_feat)]
				eigenvalues = [0.0 for i in range(n_feat)]
				eigenvectors = [[0.0 for i in range(n_feat)] for i in range(n_feat)]
				means = [0.0 for i in range(n_feat)]
				timescales = [0.0 for i in range(n_feat)]
			#
		#
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



