#!usr/bin/python

from argparse import ArgumentParser
import os
import mdtraj as md
from traj.load_features import traj_distances, traj_dihedrals, traj_dihedral_classes, load_feature_file
from tica.tica_optimize import in_equil, find_components, make_tica_opt
from output.plot_features import plot_slopes, plot_timescales, plot_density, manageable, plot_hmm, plot_estimate, plot_popl_1d, plot_estimate_1d, plot_BICs
from hmm.hmm_optimize import md_populations, make_hmm, hmm_arbitrary, make_hmm_fullauto
from hmm.inheritance import combine_hmm_1d
from checkpoint.io_results import print_tica, print_tica_projs, print_centroids
import numpy as np
import math


########################################################################
#																	   #
#					 		  Parse ARGV							   #
#																	   #
########################################################################


parser = ArgumentParser()

parser.add_argument("-i", "--input", nargs='*')
parser.add_argument("-p", "--topol")
parser.add_argument("-f", "--features") # .. # Feature file
parser.add_argument("-t", "--time") # .. # Time between frames in ps
parser.add_argument("-o", "--output")

args = vars(parser.parse_args())

trajins = args["input"]
topin = args["topol"]
f_file = args["features"]
timeskip = float(args["time"])
trajout = args["output"]

# .. # Make sure the files exist

for trajin in trajins:
#
	if not os.path.isfile(trajin):
	#
		raise Exception("Trajectory {:s} does not exist.".format(trajin))
	#
#

if not (topin is None):
#
	if not os.path.isfile(topin):
	#
		raise Exception("Topology {:s} does not exist.".format(topin))
	#
#

if not os.path.isfile(f_file):
#
	raise Exception("Feature file {:s} does not exist.".format(f_file))
#


########################################################################
#																	   #
#				  Load feature file and trajectory(ies)				   #
#																	   #
########################################################################


distances, dihedrals, dihedral_classes, exclude, keep, no_chain = load_feature_file(f_file)

max_len = 0

trajectories = []
trajectory_atoms = None

for trajin in trajins:
#
	features = None
	feature_atoms = None
	
	print("Loading trajectory {:s}...".format(trajin))
	
	traj1 = None
	
	if topin is None:
	#
		traj1 = md.load(trajin) # !! # The "topology" (the .gro file) must be included if the trajectory comes from GROMACS
	#
	else:
	#
		traj1 = md.load(trajin, top=topin) # !! # The "topology" (the .gro file) must be included if the trajectory comes from GROMACS
	#
	
	if not (keep is None):
	#
		traj1.atom_slice(traj1.top.select(keep), inplace=True)
	#
	
	if not (exclude is None):
	#
		traj1.atom_slice(traj1.top.select(exclude), inplace=True)
	#
	
	#print(instance.__class__.__name__) to print class name
	
	if not (distances is None):
	#
		calc_distances = traj_distances(traj1, distances, nochain=no_chain)
		
		if features is None:
		#
			features = calc_distances
			
			if trajectory_atoms is None:
			#
				feature_atoms = distances
			#
		#
		else:
		#
			features = np.append(features, calc_distances, axis=1)
			
			if trajectory_atoms is None:
			#
				feature_atoms = np.append(feature_atoms, distances)
			#
		#
	#
	
	if not (dihedrals is None):
	#
		calc_dihedrals = traj_dihedrals(traj1, dihedrals, nochain=no_chain)
		
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
	
	if not (dihedral_classes is None):
	#
		calc_dihedral_classes, dihedral_class_ids = traj_dihedral_classes(traj1, dihedral_classes)
		
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
	
	trajectories.append(features)
	
	if trajectory_atoms is None:
	#
		trajectory_atoms = feature_atoms
	#
#


########################################################################
#																	   #
#	 Perform a time-structure Independent Component Analysis (tICA)	   #
#	 			  Determine if equilibrium is reached				   #
#	 			 Determine optimal number of components				   #
#																	   #
########################################################################


print("Performing tICA...")

tica_tot, equil, equil_dists, n_comp, usable_comps, frameskip = make_tica_opt(trajectories, timeskip)

plot_timescales(tica_tot, usable_comps, "{:s}_timescales.pdf".format(trajout))
plot_slopes(tica_tot, usable_comps, "{:s}_slopes.pdf".format(trajout))

print_tica(tica_tot, equil, equil_dists, n_comp, usable_comps, frameskip, "{:s}_tica.dat".format(trajout))

ic_projs = tica_tot.transform(trajectories)

print_tica_projs(ic_projs, "{:s}_projs.dat".format(trajout))

all_hmms, all_BICs, all_BIC_mins, dims, pos = make_hmm_fullauto(ic_projs, equil_dists, n_comp, tica_tot.n_observations_)

hmm_opti = all_hmms[pos]

BICs = all_BICs[pos]

ic_slice = [ic_projs[0][:,:dims]]

for i in range(1, len(ic_projs)):
#
	ic_slice = [ic_projs[i][:,:dims]]
#

state_centers, mean_approx = hmm_opti.draw_centroids(ic_slice)

print("HMM with {:d} dimensions and {:d} states :".format(len(hmm_opti.means_[0]), len(hmm_opti.means_)))
print("Means :")
print(hmm_opti.means_)
print("Vars :")
print(hmm_opti.vars_)
print("Timescales :")
print(hmm_opti.timescales_)
print("Populations :")
print(hmm_opti.populations_)
print("Transition matrix :")
print(hmm_opti.transmat_)
print("Centroids :")
print(state_centers)

#plot_BICs(BICs, "{:s}_EVecs_{:d}_{:d}_BICs.pdf".format(trajout, 2*i+1, 2*i+2))

print_centroids(hmm_opti, trajins, topin, timeskip, state_centers, trajout)

for i in range(0, len(hmm_opti.means_[0]), 2):
#
	mean_slice = None
	var_slice = None
	ic_split = None
	
	if i == 0:
	#
		mean_slice = hmm_opti.means_[:,:2]
		var_slice = hmm_opti.vars_[:,:2]
		
		ic_split = [ic_projs[0][:,:2]]
		for i in range(1, len(ic_projs)):
		#
			ic_split = [ic_projs[i][:,:2]]
		#
	#
	elif i+2 == len(hmm_opti.means_[0]):
	#
		mean_slice = hmm_opti.means_[:,i-1:]
		var_slice = hmm_opti.vars_[:,i-1:]
		
		ic_split = [ic_projs[0][:,i-1:i+2]]
		for i in range(1, len(ic_projs)):
		#
			ic_split = [ic_projs[i][:,i-1:i+2]]
		#
	#
	else:
	#
		mean_slice = hmm_opti.means_[:,i-1:i+2]
		var_slice = hmm_opti.vars_[:,i-1:i+2]
		
		ic_split = [ic_projs[0][:,i-1:i+2]]
		for i in range(1, len(ic_projs)):
		#
			ic_split = [ic_projs[i][:,i-1:i+2]]
		#
	#
	
	plot_hmm(ic_split, mean_slice, var_slice, manageable(tica_tot.n_observations_), "{:s}_EVecs_{:d}_{:d}_clusters.pdf".format(trajout, i+1, i+2))
	
	plot_estimate(mean_slice, var_slice, hmm_opti.populations_, "{:s}_EVecs_{:d}_{:d}_estimate.pdf".format(trajout, i+1, i+2))
#
