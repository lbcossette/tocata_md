from msmbuilder.hmm import GaussianHMM as GHMM
import numpy as np
import math
from scipy.spatial.distance import mahalanobis as mahadist

#def popl_validate(frames1, means1, vars1, popls1):
#
	
#

def prob_element(frame, means, ivars, popls, cluster):
#
	dividend = popls[cluster] * np.exp(0.5 * mahadist(frame, means[cluster], ivars[cluster])**2) / math.sqrt(((2*math.pi)**len(means[cluster])) * np.linalg.det(ivars[cluster]))
	
	divisor = 0.0
	
	for i in range(len(popls)):
	#
		divisor += popls[i] * np.exp(0.5 * mahadist(frame, means[i], ivars[i])**2) / math.sqrt(((2*math.pi)**len(means[i])) * np.linalg.det(ivars[i]))
	#
	
	return dividend / divisor
#

def popl_intersect(frames1, means1, vars1, popls1, frames2, means2, vars2, popls2):
#
	all_trajs_1 = frames1[0]
	all_trajs_2 = frames2[0]
	
	intersect = np.zeros([len(means1), len(means2)])
	
	for i in range(len(frames1)):
	#
		all_trajs_1 = np.append(all_trajs_1, frames1[i], axis=0)
		all_trajs_2 = np.append(all_trajs_2, frames2[i], axis=0)
	#
	
	ivars1_full = np.empty([len(vars1),len(vars1[0]),len(vars1[0])])
	ivars2_full = np.empty([len(vars2),len(vars2[0]),len(vars2[0])])
	
	for i in range(len(vars1)):
	#
		ivars1_full[i] = np.diag(np.divide(np.array([1.0]*len(vars1[i])), vars1[i]))
	#
	
	for i in range(len(vars2)):
	#
		ivars2_full[i] = np.diag(np.divide(np.array([1.0]*len(vars2[i])), vars2[i]))
	#
	
	print("Inverse covariance matrices :")
	print(ivars1_full)
	print(ivars2_full)
	
	probs1 = np.empty([len(all_trajs_1),len(means1)])
	probs2 = np.empty([len(all_trajs_2),len(means2)])
	
	for i in range(len(all_trajs_1)):
	#
		for j in range(len(means1)):
		#
			#print("Frame {:d}, cluster{:d}, level 1 :".format(i,j))
			
			probs1[i][j] = prob_element(all_trajs_1[i], means1, ivars1_full, popls1, j)
			
			#print(probs1[i][j])
		#
		
		for j in range(len(means2)):
		#
			#print("Frame {:d}, cluster{:d}, level 2 :".format(i,j))
			
			probs2[i][j] = prob_element(all_trajs_2[i], means2, ivars2_full, popls2, j)
			
			#print(probs2[i][j])
		#
	#
	
	del all_trajs_1
	del all_trajs_2
	del ivars1_full
	del ivars2_full
	
	total_pop = 0.0
	
	for i in range(len(intersect)):
	#
		for j in range(len(intersect[0])):
		#
			for k in range(len(probs1)):
			#
				intersect[i][j] += probs1[k][i] * probs2[k][j]
			#
			
			intersect[i][j] /= len(probs1)
			
			total_pop += intersect[i][j]
		#
	#
	
	print("Total population : {:f}".format(total_pop))
	
	return intersect
#


def combine_hmm_1d(ic_splits, all_hmms):
#
	print("Rebuilding higher-dimensional HMMs...")
	
	rebuild_max = len(all_hmms)
	
	if rebuild_max % 2 == 1:
	#
		rebuild_max -= 1
	#
	
	means_rebuilt = [None]*int(rebuild_max/2)
	vars_rebuilt = [None]*int(rebuild_max/2)
	popls_rebuilt = [None]*int(rebuild_max/2)
	
	for i in range(len(means_rebuilt)):
	#
		print("Constructing intersection containers {:d}...".format(i+1))
		
		means_rebuilt[i] = [[None]*(2*len(all_hmms[2*i].means_[0]))]*(len(all_hmms[2*i].means_)*len(all_hmms[2*i+1].means_))
		vars_rebuilt[i] = [[None]*(2*len(all_hmms[2*i].means_[0]))]*(len(all_hmms[2*i].means_)*len(all_hmms[2*i+1].means_))
		popls_rebuilt[i] = [None]*(len(all_hmms[2*i].means_)*len(all_hmms[2*i+1].means_))
		
		print(len(all_hmms[2*i].means_)*len(all_hmms[2*i+1].means_))
		print(len(means_rebuilt[i]))
		print(len(vars_rebuilt[i]))
		print(len(popls_rebuilt[i]))
		
		print("Calculating joint populations...")
		
		popls_temp = popl_intersect(ic_splits[2*i], all_hmms[2*i].means_, all_hmms[2*i].vars_, all_hmms[2*i].populations_, ic_splits[2*i+1], all_hmms[2*i+1].means_, all_hmms[2*i+1].vars_, all_hmms[2*i+1].populations_)
		
		print(popls_temp)
		
		print("Assembly...")
		
		for j in range(len(all_hmms[2*i].means_)):
		#
			for k in range(len(all_hmms[2*i+1].means_)):
			#
				temp_means = all_hmms[2*i].means_[j]
				temp_means = np.append(temp_means, all_hmms[2*i+1].means_[k])
				
				print(j*len(all_hmms[2*i+1].means_)+k)
				
				print(temp_means)
				
				means_rebuilt[i][j*len(all_hmms[2*i+1].means_)+k] = temp_means
				
				del temp_means
				
				temp_vars = all_hmms[2*i].vars_[j]
				temp_vars = np.append(temp_vars, all_hmms[2*i+1].vars_[k])
				
				print(temp_vars)
				
				vars_rebuilt[i][j*len(all_hmms[2*i+1].means_)+k] = temp_vars
				
				del temp_vars
				
				popls_rebuilt[i][j*len(all_hmms[2*i+1].means_)+k] = popls_temp[j][k]
				
				print(popls_temp[j][k])
			#
		#
	#
	
	for i in range(len(means_rebuilt)):
	#
		print("Reconstitution {:d} :".format(i))
		
		for j in range(len(all_hmms[2*i].means_)):
		#
			for k in range(len(all_hmms[2*i+1].means_)):
			#
				print("Intersection of {:d} and {:d} :".format(j,k))
				
				print("Means:")
				print(means_rebuilt[i][j*len(all_hmms[2*i+1].means_)+k])
				
				print("Vars:")
				print(vars_rebuilt[i][j*len(all_hmms[2*i+1].means_)+k])
				
				print("Populations:")
				print(popls_rebuilt[i][j*len(all_hmms[2*i+1].means_)+k])
			#
		#
	#
	
	return means_rebuilt, vars_rebuilt, popls_rebuilt
#
