import re
import sys
import math
import numpy as np
import mdtraj as md

########################################################################
#																	   #
#			      	   	    Extract res info    					   #
#																	   #
########################################################################

def extract_res_info(topin, atom_nums, d_type):
#
	atom_infos = None
	
	for atom_num in atom_nums:
	#
		if re.match("phi", d_type):
		#
			if atom_infos is None:
			#
				atom_infos = np.array(["{:s}.{:d}.{:s}".format(topin.atom(atom_num[3]).residue.name, topin.atom(atom_num[3]).residue.resSeq, d_type)])
			#
			else:
			#
				atom_infos = np.append(atom_infos, ["{:s}.{:d}.{:s}".format(topin.atom(atom_num[3]).residue.name, topin.atom(atom_num[3]).residue.resSeq, d_type)])
			#
		#
		else:
		#
			if atom_infos is None:
			#
				atom_infos = np.array(["{:s}.{:d}.{:s}".format(topin.atom(atom_num[0]).residue.name, topin.atom(atom_num[0]).residue.resSeq, d_type)])
			#
			else:
			#
				atom_infos = np.append(atom_infos, ["{:s}.{:d}.{:s}".format(topin.atom(atom_num[0]).residue.name, topin.atom(atom_num[0]).residue.resSeq, d_type)])
			#
		#
	#
	
	return atom_infos
#

########################################################################
#																	   #
#					   		Calculate distances						   #
#																	   #
########################################################################

def traj_distances(trajin, dists, nochain=False):
#
	print("Loading specific distances...")
	
	duets = None
	
	if nochain:
	#
		for dist in dists:
		#
			match = re.match("(\w+)\.(\d+)\.(\w+) (\w+)\.(\d+)\.(\w+)$", dist)
			
			if match:
			#
				sel_str = "(resname {:s} and resSeq {:s} and name {:s}) or (resname {:s} and resSeq {:s} and name {:s})".format(match.group(1), match.group(2), match.group(3), match.group(4), match.group(5), match.group(6))
				
				new_duet = np.array([trajin.top.select(sel_str)])
				
				if duets is None:
				#
					duets = new_duet
				#
				else:
				#
					duets = np.append(duets, new_duet, axis=0)
				#
			#
			else:
			#
				raise Exception("Feature format \"{:s}\" is invalid in a no-chain distance measurement context.".format(dist))
			#
		#
	#
	else:
	#
		for dist in dists:
		#
			match = re.match("(\w+)\.(\w+)\.(\d+)\.(\w+) (\w+)\.(\w+)\.(\d+)\.(\w+)$", dist)
			
			if match:
			#
				sel_str = "(chain {:s} and resname {:s} and resSeq {:s} and name {:s}) or (chain {:s} and resname {:s} and resSeq {:s} and name {:s})".format(match.group(1), match.group(2), match.group(3), match.group(4), match.group(5), match.group(6), match.group(7), match.group(8))
				
				new_duet = np.array([trajin.top.select(sel_str)])
				
				if duets is None:
				#
					duets = new_duet
				#
				else:
				#
					duets = np.append(duets, new_duet, axis=0)
				#
			#
			else:
			#
				raise Exception("Feature format \"{:s}\" is invalid in a distance measurement context.".format(dist))
			#
		#
	#
	
	return md.compute_distances(trajin, duets)
#

########################################################################
#																	   #
#			      	   Calculate specific dihedrals 				   #
#																	   #
########################################################################

def traj_dihedrals(trajin, dihedrals, nochain=False):
#
	print("Loading specific dihedrals...")
	
	quartets = None
	
	if nochain:
	#
		for dihedral in dihedrals:
		#
			match = re.match("(\w+)\.(\d+)\.(\w+) (\w+)\.(\d+)\.(\w+) (\w+)\.(\d+)\.(\w+) (\w+)\.(\d+)\.(\w+)$", dihedral)
			
			if match:
			#
				sel_str = "(resname {:s} and resSeq {:s} and name {:s}) or (resname {:s} and resSeq {:s} and name {:s}) or (resname {:s} and resSeq {:s} and name {:s}) or (resname {:s} and resSeq {:s} and name {:s})".format(match.group(1), match.group(2), match.group(3), match.group(4), match.group(5), match.group(6), match.group(7), match.group(8), match.group(9), match.group(10), match.group(11), match.group(12))
				
				new_quartet = np.array([trajin.top.select(sel_str)])
				
				if quartets is None:
				#
					quartets = new_quartet
				#
				else:
				#
					quartets = np.append(quartets, new_quartet, axis=0)
				#
			#
			else:
			#
				raise Exception("Feature format \"{:s}\" is invalid in a no-chain dihedral measurement context.".format(dist))
			#
		#
	#
	else:
	#
		for dihedral in dihedrals:
		#
			match = re.match("(\w+)\.(\w+)\.(\d+)\.(\w+) (\w+)\.(\w+)\.(\d+)\.(\w+) (\w+)\.(\w+)\.(\d+)\.(\w+) (\w+)\.(\w+)\.(\d+)\.(\w+)$", dihedral)
			
			if match:
			#
				sel_str = "(chain {:s} and resname {:s} and resSeq {:s} and name {:s}) or (chain {:s} and resname {:s} and resSeq {:s} and name {:s}) or (chain {:s} and resname {:s} and resSeq {:s} and name {:s}) or (chain {:s} and resname {:s} and resSeq {:s} and name {:s})".format(match.group(1), match.group(2), match.group(3), match.group(4), match.group(5), match.group(6), match.group(7), match.group(8), match.group(9), match.group(10), match.group(11), match.group(12), match.group(13), match.group(14), match.group(15), match.group(16))
				
				print(sel_str)
				
				new_quartet = np.array([trajin.top.select(sel_str)])
				
				print(new_quartet)
				
				if quartets is None:
				#
					quartets = new_quartet
				#
				else:
				#
					quartets = np.append(quartets, new_quartet, axis=0)
				#
			#
			else:
			#
				raise Exception("Feature format \"{:s}\" is invalid in a dihedral measurement context.".format(dist))
			#
		#
	#
	
	return np.sin(md.compute_dihedrals(trajin, quartets))
#

########################################################################
#																	   #
#			      	   Calculate dihedral classes					   #
#																	   #
########################################################################

def traj_dihedral_classes(trajin, dihedral_classes):
#
	dihedrals_temp = None
	dihedral_ids_temp = None
	
	for dihedral_class in dihedral_classes:
	#
		print("Loading {:s} angles...".format(dihedral_class))
		
		if re.match("phi", dihedral_class):
		#
			phi_atoms, phi = md.compute_phi(trajin)
			
			if dihedrals_temp is None:
			#
				dihedrals_temp = np.cos(phi)
				
				dihedral_ids_temp = extract_res_info(trajin.top, phi_atoms, "phi")
			#
			else:
			#
				dihedrals_temp = np.append(dihedrals_temp, np.cos(phi), axis=1)
				
				dihedral_ids_temp = np.append(dihedral_ids_temp, extract_res_info(trajin.top, phi_atoms, "phi"), axis=0)
			#
		#
		
		if re.match("psi", dihedral_class):
		#
			psi_atoms, psi = md.compute_psi(trajin)
			
			if dihedrals_temp is None:
			#
				dihedrals_temp = np.sin(psi)
				
				dihedral_ids_temp = extract_res_info(trajin.top, psi_atoms, "psi")
			#
			else:
			#
				dihedrals_temp = np.append(dihedrals_temp, np.sin(psi), axis=1)
				
				dihedral_ids_temp = np.append(dihedral_ids_temp, extract_res_info(trajin.top, psi_atoms, "psi"), axis=0)
			#
		#
		
		if re.match("chi1", dihedral_class):
		#
			chi1_atoms, chi1 = md.compute_chi1(trajin)
			
			if dihedrals_temp is None:
			#
				dihedrals_temp = np.sin(chi1)
				
				dihedral_ids_temp = extract_res_info(trajin.top, chi1_atoms, "chi1")
			#
			else:
			#
				dihedrals_temp = np.append(dihedrals_temp, np.sin(chi1), axis=1)
				
				dihedral_ids_temp = np.append(dihedral_ids_temp, extract_res_info(trajin.top, chi1_atoms, "chi1"), axis=0)
			#
		#
		
		if re.match("chi2", dihedral_class):
		#
			chi2_atoms, chi2 = md.compute_chi2(trajin)
			
			if dihedrals_temp is None:
			#
				dihedrals_temp = np.sin(chi2)
				
				dihedral_ids_temp = extract_res_info(trajin.top, chi2_atoms, "chi2")
			#
			else:
			#
				dihedrals_temp = np.append(dihedrals_temp, np.sin(chi2), axis=1)
				
				dihedral_ids_temp = np.append(dihedral_ids_temp, extract_res_info(trajin.top, chi2_atoms, "chi2"), axis=0)
			#
		#
		
		if re.match("chi3", dihedral_class):
		#
			chi3_atoms, chi3 = md.compute_chi3(trajin)
			
			if dihedrals_temp is None:
			#
				dihedrals_temp = np.sin(chi3)
				
				dihedral_ids_temp = extract_res_info(trajin.top, chi3_atoms, "chi3")
			#
			else:
			#
				dihedrals_temp = np.append(dihedrals_temp, np.sin(chi3), axis=1)
				
				dihedral_ids_temp = np.append(dihedral_ids_temp, extract_res_info(trajin.top, chi3_atoms, "chi3"), axis=0)
			#
		#
		
		if re.match("chi4", dihedral_class):
		#
			chi4_atoms, chi4 = md.compute_chi4(trajin)
			
			if dihedrals_temp is None:
			#
				dihedrals_temp = np.sin(chi4)
				
				dihedral_ids_temp = extract_res_info(trajin.top, chi4_atoms, "chi4")
			#
			else:
			#
				dihedrals_temp = np.append(dihedrals_temp, np.sin(chi4), axis=1)
				
				dihedral_ids_temp = np.append(dihedral_ids_temp, extract_res_info(trajin.top, chi4_atoms, "chi4"), axis=0)
			#
		#
		
		if dihedral_class is None:
		#
			raise Exception("Feature format \"{:s}\" is invalid in a dihedral class measurement context.".format(dihedral_class))
		#
	#
	
	return (dihedrals_temp, dihedral_ids_temp)
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
