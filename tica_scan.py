#!usr/bin/python

from argparse import ArgumentParser
import os
import mdtraj as md
from analysis.compare import compare_ticas
from traj.load_features import load_trajs, load_feature_file
from traj.manip import trim_trajs, compute_features, get_sequence, get_positions, trim_make
from tica.tica_optimize import make_tica_opt
from output.plot_features import plot_scan
from checkpoint.io_results import print_feat, print_scan
import numpy as np


########################################################################
#																	   #
#					 		  Parse ARGV							   #
#																	   #
########################################################################


parser = ArgumentParser()

parser.add_argument("-i", "--input", nargs='*')
parser.add_argument("-p", "--topol")
parser.add_argument("-c", "--cover")
parser.add_argument("-f", "--features") # .. # Feature file
parser.add_argument("-t", "--time") # .. # Time between frames in ps
parser.add_argument("-o", "--output")

args = vars(parser.parse_args())

del parser

trajins = args["input"]
topin = args["topol"]
cover = int(args["cover"])
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

if cover <= 0:
#
	raise Exception("Cover must be nonzero and positive.")
#

if trajout is None:
#
	raise Exception("No specified output.")
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

print_feat(trajectory_atoms, "{:s}_features.dat".format(trajout))

sequence = get_sequence(trajectory_atoms)

tica_tot, equil, equil_dists, n_comp, usable_comps, frameskip = make_tica_opt(trajectories, timeskip)

compares = []

for i in range(n_comp):
#
	compares.append([[],[]])
#

evecs = np.asarray(tica_tot.eigenvectors_).T.tolist()[:n_comp]

centers = []

for i in range(len(sequence) - cover + 1):
#
	print("Performing trim {:d} - {:d} ...".format(sequence[i], sequence[i + cover - 1]))
	
	trim = sequence[i:i+cover]
	
	excludes = get_positions(trim, trajectory_atoms)
	
	print(excludes)
	
	temp_trajs, temp_traj_atoms = trim_make(excludes, trajectories, trajectory_atoms)
	
	print(len(temp_traj_atoms))
	print(len(temp_trajs[0][0]))
	
	if cover % 2 == 0:
	#
		centers.append(float((sequence[i + int(cover/2)] + sequence[i + int(cover/2) - 1])) / 2.0)
	#
	else:
	#
		centers.append(sequence[i + int(cover/2)])
	#
	
	print("Performing tICA...")
	
	tica_tot_c, equil_c, equil_dists_c, n_comp_c, usable_comps_c, frameskip_c = make_tica_opt(temp_trajs, timeskip)
	
	evecs_c = np.asarray(tica_tot_c.eigenvectors_).T.tolist()[:n_comp_c]
	
	coss, rtimes = compare_ticas(evecs, tica_tot.timescales_, trajectory_atoms, evecs_c, tica_tot_c.timescales_, temp_traj_atoms)
	
	for j in range(len(coss)):
	#
		compares[j][0].append(coss[j])
		compares[j][1].append(rtimes[j])
	#
#

plot_scan(compares, centers, trajout)
print_scan(compares, centers, trajout)
