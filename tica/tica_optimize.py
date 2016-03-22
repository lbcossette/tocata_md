from msmbuilder.decomposition import tICA
import math
import numpy as np


########################################################################
#																	   #
#	 			  Locate the smallest usable timescale				   #
#																	   #
########################################################################


def get_smallest_tscale(traj_tica):
#
	for i in range(len(traj_tica.timescales_)):
	#
		if traj_tica.timescales_[i] != traj_tica.timescales_[i]:
		#
			print("There are {:d} usable components.".format(i))
			
			return i
		#
	#
#


########################################################################
#																	   #
#	 Examine whether the given set of trajectories is at equilibrium   #
#																	   #
########################################################################


def in_equil(traj_tica, max_i):
#
	max_ratio = 0
	
	max_index = 0
	
	for i in range(1, max_i):
	#
		ratio = traj_tica.timescales_[i-1] / traj_tica.timescales_[i]
		
		if ratio > max_ratio:
		#
			max_ratio = ratio
			
			max_index = i
		#
	#
	
	if max_index > 1:
	#
		return (False, max_index)
	#
	elif max_index == 1:
	#
		return (True, max_index)
	#
	else:
	#
		raise Exception("Problem trying to verify equilibrium.")
	#
#


########################################################################
#																	   #
#	 			 Check the amount of relevant components			   #
#																	   #
########################################################################


def find_components(traj_tica, max_i):
#
	# .. #	Try a slope average test
	
	avg_slope = float(traj_tica.timescales_[max_i-1] - traj_tica.timescales_[0])/float(max_i-1)
	
	for i in range(2, max_i - 2):
	#
		local_slope = (traj_tica.timescales_[i+2] - traj_tica.timescales_[i-2])/8.0 + (traj_tica.timescales_[i+2] + traj_tica.timescales_[i+1] - traj_tica.timescales_[i-1] - traj_tica.timescales_[i-2])/12.0
		
		if local_slope > avg_slope and not (local_slope > 0.0):
		#
			print("There are {:d} relevant components".format(i))
			
			return i
		#
		elif local_slope > 0:
		#
			print("Problem at {:d}".format(i))
		#
	#
#


########################################################################
#																	   #
#	 			 			  Perform tICA							   #
#																	   #
########################################################################


def make_tica_opt(trajectories_t, timeskip_t):
#
	frameskip = 0
	
	if int(math.ceil(200000.0/timeskip_t)) > len(trajectories_t[0])/2:
	#
		frameskip = int(len(trajectories_t[0])/2)
		
		print("Frameskip {:f} ns ({:d} frames) :".format(float(len(trajectories_t[0])/2)*timeskip_t/1000.0, frameskip))
	#
	else:
	#
		frameskip = int(math.ceil(200000.0/timeskip_t))
		
		print("Frameskip {:d} ns ({:d} frames) :".format(200, frameskip))
	#
	
	tica_tot_t = tICA(n_components = len(trajectories_t[0][0]), lag_time = frameskip)
	
	tica_tot_t.fit(trajectories_t)
	
	usable_comps_t = get_smallest_tscale(tica_tot_t)
	
	equil_t, equil_dists_t = in_equil(tica_tot_t, usable_comps_t)
	
	n_comp_t = find_components(tica_tot_t, usable_comps_t)
	
	return tica_tot_t, equil_t, equil_dists_t, n_comp_t, usable_comps_t, frameskip
#
