import re
import sys
import math
import numpy as np
import mdtraj as md


########################################################################
#																	   #
#			      	   	     Load trajectories    					   #
#																	   #
########################################################################


def load_trajs(trajins, topin):
#
	pretrajs = []
	
	for trajin in trajins:
	#
		print("Loading {:s}".format(trajin))
		
		# .. # Load raw information from each trajectory
		
		traj1 = None
		
		if topin is None:
		#
			traj1 = md.load(trajin) # !! # The "topology" (the .gro file) must be included if the trajectory comes from GROMACS
		#
		else:
		#
			traj1 = md.load(trajin, top=topin)
		#
		
		pretrajs.append(traj1)
	#
	
	return pretrajs
#


########################################################################
#																	   #
#			      	   	    Load feature file    					   #
#																	   #
########################################################################


def load_feature_file(f_file):
#
	feat_file = open(f_file, 'r')

	distances = None
	dihedrals = None
	dihedral_classes = None
	exclude = None
	keep = None
	
	is_distance = False
	is_dihedral = False
	is_dihedral_class = False
	is_exclude = False
	is_keep = False
	no_chain = False
	
	for line in feat_file:
	#
		print(line.strip('\n'))
	
		if re.match("Distances", line):
		#
			is_distance = True
			is_dihedral = False
			is_dihedral_class = False
			is_exclude = False
			is_keep = False
		
			continue
		#
	
		if re.match("Dihedrals", line):
		#
			is_distance = False
			is_dihedral = True
			is_dihedral_class = False
			is_exclude = False
			is_keep = False
		
			continue
		#
	
		if re.match("Dihedral classes", line):
		#
			is_distance = False
			is_dihedral = False
			is_dihedral_class = True
			is_exclude = False
			is_keep = False
		
			continue
		#
		
		if re.match("Exclude", line):
		#
			is_distance = False
			is_dihedral = False
			is_dihedral_class = False
			is_exclude = True
			is_keep = False
		
			continue
		#
		
		if re.match("Keep", line):
		#
			is_distance = False
			is_dihedral = False
			is_dihedral_class = False
			is_exclude = False
			is_keep = True
		
			continue
		#
		
		if re.match("No chain", line):
		#
			no_chain = True
		
			continue
		#
	
		if is_distance:
		#
			if distances is None:
			#
				distances = np.array([line.strip('\n')])
			#
			else:
			#
				distances = np.append(distances, [line.strip('\n')])
			#
		
			continue
		#
	
		if is_dihedral:
		#
			if dihedrals is None:
			#
				dihedrals = np.array([line.strip('\n')])
			#
			else:
			#
				dihedrals = np.append(dihedrals, [line.strip('\n')])
			#
		
			continue
		#
	
		if is_dihedral_class:
		#
			if dihedral_classes is None:
			#
				dihedral_classes = np.array([line.strip('\n')])
			#
			else:
			#
				dihedral_classes = np.append(dihedral_classes, [line.strip('\n')])
			#
		
			continue
		#
		
		if is_exclude:
		#
			exclude = "not ({:s})".format(line.strip('\n'))
		#
		
		if is_keep:
		#
			keep = line.strip('\n')
		#
	#
	
	return (distances, dihedrals, dihedral_classes, exclude, keep, no_chain)
#
