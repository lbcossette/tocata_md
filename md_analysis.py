#!usr/bin/python

from argparse import ArgumentParser
import os
import mdtraj as md
from traj.load_features import load_trajs, load_feature_file
from traj.manip import trim_trajs, compute_features, get_sequence
from tica.tica_optimize import make_tica_opt
from output.plot_features import plot_slopes, plot_timescales, manageable, plot_hmm, plot_estimate, plot_BICs, plot_trajs
from hmm.hmm_optimize import hmm_arbitrary, make_hmm_fullauto
from analysis.inheritance import inheritance
from checkpoint.io_results import print_tica, print_tica_projs, print_centroids, print_hmm, print_feat
import numpy as np


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
parser.add_argument("-s", "--structure")
parser.add_argument("-a", "--arbitrary", nargs='*') # .. # For arbitrary hmm construction (n_dimensions, n_clusters)
parser.add_argument("-r", "--resume", nargs='*') # .. # To resume hmm construction or override exclusion criteria (n_dimensions, n_clusters)
parser.add_argument("-d", "--data")

args = vars(parser.parse_args())

del parser

trajins = args["input"]
topin = args["topol"]
f_file = args["features"]
timeskip = float(args["time"])
trajout = args["output"]
structin = args["structure"]
arbitrary = args["arbitrary"]
resume = args["resume"]
print_data = args["data"]

toprint = True

if not (print_data is None):
#
	if print_data == "0":
	#
		toprint = False
	#
	elif not (print_data == "1"):
	#
		raise Exception("Could not parse -d argument.")
	#
#

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
	
	if not os.path.isfile(structin):
	#
		raise Exception("Structure {:s} does not exist.".format(structin))
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


# .. # Load feature file for features to extract

distances, dihedrals, dihedral_classes, exclude, keep, no_chain = load_feature_file(f_file)

del f_file

# .. # Load the selected features from the trajectory

pretrajs = load_trajs(trajins, topin)

pretrajs = trim_trajs(pretrajs, exclude, keep)

trajectories, trajectory_atoms = compute_features(pretrajs, distances, dihedrals, dihedral_classes, no_chain)

sequence = get_sequence(trajectory_atoms)

#if toprint:
#
	#plot_trajs(trajectories, trajectory_atoms, sequence, trajins, timeskip, trajout)
#

print_feat(trajectory_atoms, "{:s}_features.dat".format(trajout))

########################################################################
#																	   #
#	 Perform a time-structure Independent Component Analysis (tICA)	   #
#	 			  Determine if equilibrium is reached				   #
#	 			 Determine optimal number of components				   #
#																	   #
########################################################################

print(np.shape(trajectories[0]))

print("Performing tICA...")

tica_tot, equil, equil_dists, n_comp, usable_comps, frameskip = make_tica_opt(trajectories, timeskip)

# .. # Plot the timescales of the tICA and their estimated derivatives, then print all relevant tICA information

plot_timescales(tica_tot, usable_comps, "{:s}_timescales.pdf".format(trajout))
plot_slopes(tica_tot, usable_comps, "{:s}_slopes.pdf".format(trajout))

if toprint:
#
	print_tica(tica_tot, equil, equil_dists, n_comp, usable_comps, frameskip, "{:s}_tica.dat".format(trajout))
#

# .. # Project all trajectory data on independent components, then print the projections

ic_projs = tica_tot.transform(trajectories)

n_obs = tica_tot.n_observations_

del trajectories
del trajectory_atoms
del tica_tot
del equil
del usable_comps

if toprint:
#
	print_tica_projs(ic_projs, n_comp, "{:s}_projs.dat".format(trajout))
#

hmm_opti = None
ic_slice = None

if arbitrary is None:
#
	all_hmms, all_BICs, all_BIC_mins, dims, pos, all_dims = make_hmm_fullauto(ic_projs, equil_dists, n_comp, n_obs, resume)
	
	hmm_opti = all_hmms[pos]
	
	BICs = all_BICs[pos]
	
	ic_slice = [ic_projs[0][:,:dims]]
	for i in range(1, len(ic_projs)):
	#
		ic_slice.append(ic_projs[i][:,:dims])
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
	
	for i in range(len(hmm_opti.transmat_)):
	#
		print(", ".join(map(str, hmm_opti.transmat_[i])))
	#
	
	print("Centroids :")
	print(state_centers)
	
	if toprint:
	#
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
				for j in range(1, len(ic_projs)):
				#
					ic_split.append(ic_projs[j][:,:2])
				#
			#
			elif i+2 == len(hmm_opti.means_[0]):
			#
				mean_slice = hmm_opti.means_[:,i:]
				var_slice = hmm_opti.vars_[:,i:]
				
				ic_split = [ic_projs[0][:,i:i+2]]
				for j in range(1, len(ic_projs)):
				#
					ic_split.append(ic_projs[j][:,i:i+2])
				#
			#
			else:
			#
				mean_slice = hmm_opti.means_[:,i:i+2]
				var_slice = hmm_opti.vars_[:,i:i+2]
				
				ic_split = [ic_projs[0][:,i:i+2]]
				for j in range(1, len(ic_projs)):
				#
					ic_split.append(ic_projs[j][:,i:i+2])
				#
			#
			
			plot_hmm(ic_split, mean_slice, var_slice, manageable(n_obs), "{:s}_EVecs_{:d}_{:d}_clusters.pdf".format(trajout, i+1, i+2))
			
			plot_estimate(mean_slice, var_slice, hmm_opti.populations_, "{:s}_EVecs_{:d}_{:d}_estimate.pdf".format(trajout, i+1, i+2))
		#
			
		print_centroids(hmm_opti, trajins, structin, timeskip, state_centers, trajout)
		
		#plot_BICs(all_BICs, all_dims, "{:s}_BICs.pdf".format(trajout, 2*i+1, 2*i+2))
	#
#
else:
#
	ic_slice = [ic_projs[0][:,:int(arbitrary[0])]]
	for i in range(1, len(ic_projs)):
	#
		ic_slice.append(ic_projs[i][:,:int(arbitrary[0])])
	#
	
	hmm_opti, BIC = hmm_arbitrary(ic_slice, int(arbitrary[0]), int(arbitrary[1]), n_obs)
	
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
	
	for i in range(len(hmm_opti.transmat_)):
	#
		print(", ".join(map(str, hmm_opti.transmat_[i])))
	#
	
	print("Centroids :")
	print(state_centers)
	print("BIC :")
	print(BIC)
	
	if toprint:
	#
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
				for j in range(1, len(ic_projs)):
				#
					ic_split.append(ic_projs[j][:,:2])
				#
			#
			elif i+2 == len(hmm_opti.means_[0]):
			#
				mean_slice = hmm_opti.means_[:,i:]
				var_slice = hmm_opti.vars_[:,i:]
				
				ic_split = [ic_projs[0][:,i:i+2]]
				for j in range(1, len(ic_projs)):
				#
					ic_split.append(ic_projs[j][:,i:i+2])
				#
			#
			else:
			#
				mean_slice = hmm_opti.means_[:,i:i+2]
				var_slice = hmm_opti.vars_[:,i:i+2]
				
				ic_split = [ic_projs[0][:,i:i+2]]
				for j in range(1, len(ic_projs)):
				#
					ic_split.append(ic_projs[j][:,i:i+2])
				#
			#
			
			plot_hmm(ic_split, mean_slice, var_slice, manageable(n_obs), "{:s}_EVecs_{:d}_{:d}_clusters.pdf".format(trajout, i+1, i+2))
			
			plot_estimate(mean_slice, var_slice, hmm_opti.populations_, "{:s}_EVecs_{:d}_{:d}_estimate.pdf".format(trajout, i+1, i+2))
		#
		
		print_centroids(hmm_opti, trajins, structin, timeskip, state_centers, trajout)
	#
#

print_hmm(hmm_opti, "{:s}_hmm.dat".format(trajout))
