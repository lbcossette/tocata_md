import re
import sys
import math
import numpy as np
import mdtraj as md
from copy import deepcopy


########################################################################
#																	   #
#			      	   	    Extract res info    					   #
#																	   #
########################################################################

def extract_res_info(topin, atom_nums, d_type): # <> # Currenty only works on "no chain" mode. To be fixed eventually.
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
#			      	   	Get positions to exclude   					   #
#																	   #
########################################################################


def get_positions(trim, trajectory_atoms):
#
	excludes = []
	temp_traj_atoms = []
	
	for i in range(len(trim)):
	#
		for j in range(len(trajectory_atoms)):
		#
			if int(re.match("\w+\.(\d+)", trajectory_atoms[j]).group(1)) == trim[i]:
			#
				excludes.append(j)
			#
		#
	#
	
	ord_excludes = np.sort(excludes)
	
	return ord_excludes
#


########################################################################
#																	   #
#			  Trim trajectories and trajectory atoms list			   #
#																	   #
########################################################################


def trim_make(excludes, trajectories, trajectory_atoms):
#
	new_trajs = []
	
	new_traj_atoms = trajectory_atoms[:excludes[0]]
	
	for i in range(1,len(excludes)):
	#
		new_traj_atoms = np.append(new_traj_atoms, trajectory_atoms[excludes[i-1]+1:excludes[i]])
	#
	
	new_traj_atoms = np.append(new_traj_atoms, trajectory_atoms[excludes[-1]+1:])
	
	for trajectory in trajectories:
	#
		new_traj = trajectory[...,:excludes[0]]
		
		for i in range(1, len(excludes)):
		#
			new_traj = np.append(new_traj, trajectory[...,excludes[i-1]+1:excludes[i]], axis=1)
		#
		
		new_traj = np.append(new_traj, trajectory[...,excludes[-1]+1:], axis=1)
		
		new_trajs.append(new_traj)
	#
	
	return new_trajs, new_traj_atoms
#


########################################################################
#																	   #
#			      	   	       Get sequence	    					   #
#																	   #
########################################################################


def get_sequence(features): # <> # Currently only works on "no chain" mode. To be fixed eventually.
#
	sequence = []
	
	for feature in features:
	#
		m = re.match("\w+\.(\d+)", feature)
		
		pos = int(m.group(1))
		
		if not sequence:
		#
			sequence.append(pos)
		#
		elif pos > sequence[-1]:
		#
			sequence.append(pos)
		#
		elif pos < sequence[0]:
		#
			new_seq = [pos]
			
			for ele in sequence:
			#
				new_seq.append(ele)
			#
			
			sequence = new_seq
		#
		else:
		#
			for i in range(1, len(sequence)):
			#
				if pos > sequence[i-1] and pos < sequence[i]:
				#
					seq_less = deepcopy(sequence[:i])
					
					seq_more = deepcopy(sequence[i:])
					
					sequence = []
					
					for s_less in seq_less:
					#
						sequence.append(s_less)
					#
					
					sequence.append(pos)
					
					for s_more in seq_more:
					#
						sequence.append(s_more)
					#
				#
			#
		#
	#
	
	return sequence
#


########################################################################
#																	   #
#					   		Calculate distances						   #
#																	   #
########################################################################


def traj_distances(trajin, dists, nochain=False):
#
	#print("Loading specific distances...")
	
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
	#print("Loading specific dihedrals...")
	
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
		#print("Loading {:s} angles...".format(dihedral_class))
		
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
#			      	   	    Trim trajectories    					   #
#																	   #
########################################################################


def trim_trajs(pretrajs, exclude, keep):
#
	print("Trimming trajectories...")
	
	for pretraj in pretrajs:
	#
		#print(pretraj)
		
		if not (keep is None):
		#
			pretraj.atom_slice(pretraj.top.select(keep), inplace=True)
		#
		
		if not (exclude is None):
		#
			pretraj.atom_slice(pretraj.top.select(exclude), inplace=True)
		#
		
		#print(pretraj)
	#
	
	return pretrajs
#


########################################################################
#																	   #
#			      	   	 Trim/copy trajectories    					   #
#																	   #
########################################################################


def trim_trajs_copy(pretrajs, exclude, keep):
#
	print("Trimming trajectories...")
	
	newtrajs = []
	
	for pretraj in pretrajs:
	#
		#print(pretraj)
		
		newtraj = None
		
		if not (keep is None):
		#
			newtraj = pretraj.atom_slice(pretraj.top.select(keep), inplace=False)
		#
		
		if not (exclude is None):
		#
			if not (newtraj is None):
			#
				newtraj.atom_slice(newtraj.top.select(exclude), inplace=True)
			#
			else:
			#
				newtraj = pretraj.atom_slice(pretraj.top.select(exclude), inplace=False)
			#
		#
		
		#print(pretraj)
		#print(newtraj)
		
		newtrajs.append(newtraj)
	#
	
	return newtrajs
#


########################################################################
#																	   #
#			      	   	     Compute features    					   #
#																	   #
########################################################################


def compute_features(pretrajs, distances, dihedrals, dihedral_classes, no_chain):
#
	print("Computing features...")
	
	trajectories = []
	trajectory_atoms = None
	
	for pretraj in pretrajs:
	#
		print(pretraj)
		
		features = None
		feature_atoms = None
		
		# .. # Extract specific distances from trajectory
		
		if not (distances is None):
		#
			calc_distances = traj_distances(pretraj, distances, nochain=no_chain)
			
			features = calc_distances
			
			if trajectory_atoms is None:
			#
				feature_atoms = distances
			#
		#
		
		# .. # Extract specific dihedrals from trajectory
		
		if not (dihedrals is None):
		#
			calc_dihedrals = traj_dihedrals(pretraj, dihedrals, nochain=no_chain)
			
			if features is None:
			#
				features = calc_dihedrals
				
				if trajectory_atoms is None:
				#
					feature_atoms = dihedrals
				#
			#
			else:
			#
				features = np.append(features, calc_dihedrals, axis=1)
				
				if trajectory_atoms is None:
				#
					feature_atoms = np.append(feature_atoms, dihedrals)
				#
			#
		#
		
		# .. # Extract all specified dihedral classes (phi, psi, chi1, chi2, chi3, chi4) from trajectory
		
		if not (dihedral_classes is None):
		#
			calc_dihedral_classes, dihedral_class_ids = traj_dihedral_classes(pretraj, dihedral_classes)
			
			if features is None:
			#
				features = calc_dihedral_classes
				
				if trajectory_atoms is None:
				#
					feature_atoms = dihedral_class_ids
				#
			#
			else:
			#
				features = np.append(features, calc_dihedral_classes, axis=1)
				
				if trajectory_atoms is None:
				#
					feature_atoms = np.append(feature_atoms, dihedral_class_ids)
				#
			#
		#
		
		print("Done.")
		
		# .. # Append the transformed trajectory to the existing list of trajectories
		
		trajectories.append(features)
		
		if trajectory_atoms is None:
		#
			trajectory_atoms = feature_atoms
		#
	#
	
	return trajectories, trajectory_atoms
#
