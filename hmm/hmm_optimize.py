from msmbuilder.hmm import GaussianHMM as GHMM
import numpy as np
from copy import deepcopy
from scipy.stats import chi2
import math


########################################################################
#																	   #
#	 	   Calculate the actual populations of the trajectory		   #
#																	   #
########################################################################


def md_populations(ic_slice, dimension):
#	
	print(len(ic_slice[0][0]))
	
	if len(ic_slice[0][0]) == 1:
	#
		print("Switching to 1D...")
		
		return md_populations_1d(ic_slice, dimension)
	#
	
	print("Calculating populations...")
	
	all_trajs = ic_slice[0]
	
	for i in range(1, len(ic_slice)):
	#
		all_trajs = np.append(all_trajs, ic_slice[i], axis=0)
	#
	
	x_range = [all_trajs[0][0], all_trajs[0][0]]
	
	y_range = [all_trajs[0][1], all_trajs[0][1]]
	
	for i in range(len(all_trajs)):
	#
		if all_trajs[i][1] > y_range[1]:
		#
			y_range[1] = all_trajs[i][1]
		#
		elif all_trajs[i][1] < y_range[0]:
		#
			y_range[0] = all_trajs[i][1]
		#
		
		if all_trajs[i][0] > x_range[1]:
		#
			x_range[1] = all_trajs[i][0]
		#
		elif all_trajs[i][0] < x_range[0]:
		#
			x_range[0] = all_trajs[i][0]
		#
	#
	
	populations_md = np.array([[0.0]*dimension]*dimension)
	
	for i in range(len(all_trajs)):
	#
		x_pos = int((all_trajs[i][0] - x_range[0])/(x_range[1] - x_range[0]) * float(dimension))
		
		y_pos = int((all_trajs[i][1] - y_range[0])/(y_range[1] - y_range[0]) * float(dimension))
		
		if x_pos != dimension and y_pos != dimension:
		#
			populations_md[x_pos][y_pos] += 1.0
		#
		elif x_pos == dimension and y_pos != dimension:
		#
			populations_md[dimension - 1][y_pos] += 1.0
		#
		elif x_pos != dimension and y_pos == dimension:
		#
			populations_md[x_pos][dimension - 1] += 1.0
		#
		else:
		#
			populations_md[dimension - 1][dimension - 1] += 1.0
		#
	#
	
	for i in range(dimension):
	#
		for j in range(dimension):
		#
			populations_md[i][j] /= float(len(all_trajs))
		#
	#
	
	print("Done.")
	
	return (populations_md, [x_range, y_range], 0.5 / float(len(all_trajs)), float(len(all_trajs)))
#


########################################################################
#																	   #
#	 		  Calculate the estimated populations of HMM			   #
#																	   #
########################################################################


def hmm_populations(slice_hmm, dimension, ranges, min_pop):
#
	if len(slice_hmm.means_[0]) == 1:
	#
		print("Switching to 1D...")
		
		return hmm_populations_1d(slice_hmm, dimension, ranges, min_pop)
	#
	
	print("Calculating estimated populations...")
	
	# <> # Due to constraints related to the "nquad" function and the usage of hmm_dist, it is better to include it in hmm_populations
	# <> # This constraint is not true anymore, but for now the function is left "as is". It is not useful to make it independent yet.
	
	def hmm_dist(x, y):
	#
		value = 0.0
		
		for i in range(len(slice_hmm.means_)):
		#
			arg_exp = 0.5*((x - slice_hmm.means_[i][0])**2 / slice_hmm.vars_[i][0] + (y - slice_hmm.means_[i][1])**2 / slice_hmm.vars_[i][1])
			
			divisor = 2 * np.pi * np.sqrt(slice_hmm.vars_[i][0] * slice_hmm.vars_[i][1])
			
			value += slice_hmm.populations_[i] * np.exp(-1.0 * arg_exp) / divisor
		#
		
		return value
	#
	
	populations_hmm = np.array([[None]*dimension]*dimension)
	
	total_est = 0.0
	
	for i in range(dimension):
	#
		for j in range(dimension):
		#
			x1 = ranges[0][0] + i*(ranges[0][1] - ranges[0][0])/dimension
			x2 = ranges[0][0] + (i+1)*(ranges[0][1] - ranges[0][0])/dimension
			y1 = ranges[1][0] + j*(ranges[1][1] - ranges[1][0])/dimension
			y2 = ranges[1][0] + (j+1)*(ranges[1][1] - ranges[1][0])/dimension
			
			dist_avg = hmm_dist(x1, y1) + hmm_dist(x2, y1) + hmm_dist(x1, y2) + hmm_dist(x2, y2)
			
			dist_avg /= 4.0
			
			populations_hmm[i][j] = (ranges[0][1] - ranges[0][0])/dimension * (ranges[1][1] - ranges[1][0])/dimension * dist_avg
			
			if populations_hmm[i][j] < min_pop:
			#
				populations_hmm[i][j] = min_pop
			#
			
			total_est += populations_hmm[i][j]
		#
	#
	
	for i in range(dimension):
	#
		for j in range(dimension):
		#
			populations_hmm[i][j] /= total_est
		#
	#
	
	print("Done")
	
	return populations_hmm
#


########################################################################
#																	   #
#		  	Calculate the actual populations of the trajectory		   #
#							1D test function						   #
#																	   #
########################################################################


def md_populations_1d(ic_slice, dimension):
#
	print("Calculating populations...")
	
	all_trajs = ic_slice[0]
	
	for i in range(1, len(ic_slice)):
	#
		all_trajs = np.append(all_trajs, ic_slice[i], axis=0)
	#
	
	x_range = [all_trajs[0][0], all_trajs[0][0]]
	
	for i in range(len(all_trajs)):
	#
		if all_trajs[i][0] > x_range[1]:
		#
			x_range[1] = all_trajs[i][0]
		#
		elif all_trajs[i][0] < x_range[0]:
		#
			x_range[0] = all_trajs[i][0]
		#
	#
	
	populations_md = np.array([0.0]*dimension)
	
	for i in range(len(all_trajs)):
	#
		x_pos = int((all_trajs[i][0] - x_range[0])/(x_range[1] - x_range[0]) * float(dimension))
		
		if x_pos == dimension:
		#
			populations_md[x_pos-1] += 1.0
		#
		else:
		#
			populations_md[x_pos] += 1.0
		#
	#
	
	for i in range(dimension):
	#
		populations_md[i] /= float(len(all_trajs))
	#
	
	print("Done.")
	
	return (populations_md, [x_range], 0.5 / float(len(all_trajs)), float(len(all_trajs)))
#


########################################################################
#																	   #
#	 		  Calculate the estimated populations of HMM			   #
#							1D test function						   #
#																	   #
########################################################################


def hmm_populations_1d(slice_hmm, dimension, ranges, min_pop):
#
	print("Calculating estimated populations...")
	
	def hmm_dist_1d(x):
	#
		value = 0.0
		
		for i in range(len(slice_hmm.means_)):
		#
			arg_exp = 0.5*((x - slice_hmm.means_[i][0])**2 / slice_hmm.vars_[i][0])
			
			divisor = 2 * np.pi * np.sqrt(slice_hmm.vars_[i][0])
			
			value += slice_hmm.populations_[i] * np.exp(-1.0 * arg_exp) / divisor
		#
		
		return value
	#
	
	populations_hmm = np.array([None]*dimension)
	
	total_est = 0.0
	
	for i in range(dimension):
	#
		x1 = ranges[0][0] + i*(ranges[0][1] - ranges[0][0])/dimension
		x2 = ranges[0][0] + (i+1)*(ranges[0][1] - ranges[0][0])/dimension
		
		dist_avg = hmm_dist_1d(x1) + hmm_dist_1d(x2)
		
		dist_avg /= 2.0
		
		populations_hmm[i] = (ranges[0][1] - ranges[0][0])/dimension * dist_avg
		
		if populations_hmm[i] < min_pop:
		#
			populations_hmm[i] = min_pop
		#
		
		total_est += populations_hmm[i]
	#
	
	for i in range(dimension):
	#
			populations_hmm[i] /= total_est
	#
	
	print("Done")
	
	return populations_hmm
#


########################################################################
#																	   #
#     Calculate the ratio between the Kullback-Leibler divergence 	   #
#                 and the entropy of the posterior PDF				   #
#																	   #
########################################################################


# <> # KLD and H explode on low-data limit. Outliers (slots with 1 or 2 points) are removed with the 2.97*min_pop criterion.


def hmm_KLD_H(populations_hmm, populations, min_pop):
#
	if isinstance(populations[0], float):
	#
		print("Switching to 1D...")
		
		return hmm_KLD_H_1d(populations_hmm, populations, min_pop)
	#
	
	print("Calculating Kullback-Leibler Divergence ratio...")
	
	KLD = 0.0
	H = 0.0
	
	for i in range(len(populations)):
	#
		for j in range(len(populations[0])):
		#
			if populations[i][j] >= 2.97*min_pop:
			#
				KLD += populations[i][j] * math.log(populations[i][j] / populations_hmm[i][j])
				H -= populations[i][j] * math.log(populations[i][j])
			#
		#
	#
	
	print("Kullback-Leibler Divergence :")
	print(KLD)
	print("Entropy :")
	print(H)
	
	print("Done.")
	
	return KLD / H
#


########################################################################
#																	   #
#     Calculate the ratio between the Kullback-Leibler divergence 	   #
#                 and the entropy of the posterior PDF				   #
#							1D test function						   #
#																	   #
########################################################################


# <> # KLD and H explode on low-data limit. Outliers (slots with 1 or 2 points) are removed with the 2.97*min_pop criterion.


def hmm_KLD_H_1d(populations_hmm, populations, min_pop):
#
	print("Calculating Kullback-Leibler Divergence ratio...")
	
	KLD = 0.0
	H = 0.0
	
	for i in range(len(populations)):
	#
		if populations[i] >= 0.99*min_pop:
		#
			KLD += populations[i] * math.log(populations[i] / populations_hmm[i])
			H -= populations[i] * math.log(populations[i])
		#
	#
	
	print("Kullback-Leibler Divergence :")
	print(KLD)
	print("Entropy :")
	print(H)
	
	print("Done.")
	
	return KLD / H
#


########################################################################
#																	   #
#	 		  			Build and optimize HMM						   #
#																	   #
########################################################################


def make_hmm(ic_slice, dimension, dims):
#
	#bad_clustering = True
	#while bad_clustering:
		#n_clusters += 1
	
	populations, ranges, min_pop, tot_pop = md_populations(ic_slice, dimension)
	
	n_clusters = 0
	
	best_best_hmm = None
	#best_best_KLD_H = None
	best_best_hmm_popl = None
	
	BICs = np.zeros(1)
	
	BIC_min = 0
	
	strike_max = 0
	strike_min = 0
	
	while True:
	#
		n_clusters += 1
		
		print("Trying {:d} clusters...".format(n_clusters))
		
		hmm = GHMM(n_states=n_clusters, reversible_type='transpose', thresh=1e-4, init_algo='GMM', timing=True)
		
		best_hmm = None
		#best_KLD_H = None
		best_BIC = None
		
		for tries in range(1):
		#
			print("Run {:d}...".format(tries + 1))
			
			hmm.fit(ic_slice)
			
			populations_hmm = hmm_populations(hmm, dimension, ranges, min_pop)
			
			if tot_pop > 10000.0:
			#
				BIC_hmm = (1.0 + 2.0 * float(dims)) * float(n_clusters)**2 * math.log(tot_pop) - 2.0 * hmm.fit_logprob_[0]
			#
			else:
			#
				BIC_hmm = (1.0 + 2.0 * float(dims)) * float(n_clusters) * math.log(tot_pop) - 2.0 * hmm.fit_logprob_[0]
			#
			
			#KLD_H_hmm = hmm_KLD_H(populations_hmm, populations, min_pop)
			
			if best_hmm is None:
			#
				print("Better fit hmm found with {:d} clusters".format(n_clusters))
				#print(KLD_H_hmm)
				#print(best_KLD_H)
				print("BIC :")
				print(BIC_hmm)
				print("Penalty :")
				print((1.0 + 2.0 * float(dims)) * float(n_clusters) * math.log(tot_pop))
				print((1.0 + 2.0 * float(dims)) * float(n_clusters)**2 * math.log(tot_pop))
				print((1.0 + 2.0 * float(dims)) * float(n_clusters)**3 * math.log(tot_pop))
				
				best_hmm = deepcopy(hmm)
				#best_KLD_H = KLD_H_hmm
				best_hmm_popl = populations_hmm
				best_BIC = BIC_hmm
			#
			elif BIC_hmm < best_BIC:
			#
				print("Better fit hmm found with {:d} clusters".format(n_clusters))
				#print(KLD_H_hmm)
				#print(best_KLD_H)
				print("BIC :")
				print(BIC_hmm)
				print("Penalty :")
				print((1.0 + 2.0 * float(dims)) * float(n_clusters) * math.log(tot_pop))
				print((1.0 + 2.0 * float(dims)) * float(n_clusters)**2 * math.log(tot_pop))
				print((1.0 + 2.0 * float(dims)) * float(n_clusters)**3 * math.log(tot_pop))
				print("Previous minimum :")
				print(best_BIC)
				
				best_hmm = deepcopy(hmm)
				#best_KLD_H = KLD_H_hmm
				best_hmm_popl = populations_hmm
				best_BIC = BIC_hmm
			#
		#
		
		if BICs[0] == 0.0:
		#
			BICs[0] = best_BIC
		#
		else:
		#
			BICs = np.append(BICs, [best_BIC])
		#
		
		if len(BICs) == 1:
		#
			best_best_hmm = deepcopy(best_hmm)
			best_best_hmm_popl = best_hmm_popl
		#
		elif len(BICs) > 1:
		#
			if BICs[len(BICs)-1] >= BICs[len(BICs)-2]:
			#
				print("Strike max")
			
				strike_max += 1
			#
			else:
			#
				print("Reset max")
			
				strike_max = 0
			
				if BICs[len(BICs)-1] >= BICs[BIC_min]:
				#
					if len(BICs) > 2 and BICs[len(BICs)-2] <= BICs[len(BICs)-3]:
					#
						print("Carry min")
					#
					else:
					#
						print("Strike min")
					
						strike_min += 1
					#
				#
				else:
				#
					print("Reset min")
				
					strike_min = 0
				
					BIC_min = len(BICs)-1
					
					best_best_hmm = deepcopy(best_hmm)
					best_best_hmm_popl = best_hmm_popl
				#
			#
		#
		
		if n_clusters >= 10 and (strike_max >= 3 or strike_min >= 3):
		#
			print("Reached expected true minimum")
			
			break
		#
		
		#print("Kullback-Leibler Divergence to posterior entropy ratio :")
		#print(best_KLD_H)
		
		#if best_KLD_H <= 0.005:
		#
			#print("Information loss is inferior to 0.75%")
			
			#best_best_KLD_H = best_KLD_H
			#best_best_hmm = deepcopy(best_hmm)
			#best_best_hmm_popl = best_hmm_popl
			
			#break
		#
		
		print("\n")
	#
	
	return best_best_hmm, ranges, populations, best_best_hmm_popl, BICs
#


########################################################################
#																	   #
#	 		  			Generate arbitrary HMM						   #
#																	   #
########################################################################


def hmm_arbitrary(ic_slice, n_dims, n_clusters, n_obs):
#
	BIC_hmm = 0.0
	
	print("Trying {:d} clusters...".format(n_clusters))
	
	hmm = GHMM(n_states=n_clusters, reversible_type='transpose', thresh=1e-4, init_algo='GMM', timing=True)
	
	print("Running fit...")
	
	hmm.fit(ic_slice)
	
	if n_obs > 10000:
	#
		BIC_hmm = (1.0 + 2.0 * float(n_dims)) * float(n_clusters)**2 * math.log(n_obs) - 2.0 * hmm.fit_logprob_[0]
	#
	else:
	#
		BIC_hmm = (1.0 + 2.0 * float(n_dims)) * float(n_clusters) * math.log(n_obs) - 2.0 * hmm.fit_logprob_[0]
	#
	
	return hmm, BIC_hmm
#


########################################################################
#																	   #
#	 		  		 Fully automatic HMM generation					   #
#																	   #
########################################################################


def make_hmm_fullauto(ic_projs, equil_dists, n_comp, n_obs, resume):
#
	all_hmms = [None]
	all_BICs = [None]
	all_BIC_mins = [None]
	min_dim = 0
	pos = 0
	
	if equil_dists > n_comp:
	#
		temp = equil_dists
		
		equil_dists = n_comp
		
		n_comp = equil_dists
		
		del temp
	#
	
	n_dims = equil_dists - 1
	
	while n_dims <= n_comp:
	#
		n_dims += 1
		
		n_clusters = 0
		
		if not (resume is None):
		#
			print("Resuming with {:d} dimensions and {:d} clusters.".format(int(resume[0]), int(resume[1])))
			
			n_dims = int(resume[0])
			n_clusters = int(resume[1]) - 1
			resume = None
		#
		
		print("Trying with {:d} dimensions...\n".format(n_dims))
		
		ic_slice = [ic_projs[0][:,:n_dims]]
		
		for j in range(1, len(ic_projs)):
		#
			ic_slice.append(ic_projs[j][:,:n_dims])
		#
		
		best_hmm = None
		
		BICs = np.zeros(1)
		BIC_min = 0
		
		strike_max = 0
		strike_min = 0
		
		while True:
		#
			n_clusters += 1
			
			print("Trying {:d} clusters...".format(n_clusters))
			
			hmm = GHMM(n_states=n_clusters, reversible_type='transpose', thresh=1e-4, init_algo='GMM', timing=True)
			
			print("Running fit...")
			
			hmm.fit(ic_slice)
			
			if n_obs > 10000:
			#
				BIC_hmm = (1.0 + 2.0 * float(n_dims)) * float(n_clusters)**2 * math.log(n_obs) - 2.0 * hmm.fit_logprob_[-1]
			#
			else:
			#
				BIC_hmm = (1.0 + 2.0 * float(n_dims)) * float(n_clusters) * math.log(n_obs) - 2.0 * hmm.fit_logprob_[-1]
			#
			
			if BICs[0] == 0.0:
			#
				BICs[0] = BIC_hmm
			#
			else:
			#
				BICs = np.append(BICs, [BIC_hmm])
			#
			
			if len(BICs) == 1:
			#
				best_hmm = deepcopy(hmm)
			#
			elif len(BICs) > 1:
			#
				if BICs[len(BICs)-1] >= BICs[len(BICs)-2]:
				#
					print("Strike max")
					print("Current BIC:")
					print(BIC_hmm)
					print("Current min:")
					print(BICs[BIC_min])
					
					strike_max += 1
				#
				else:
				#
					print("Reset max")
					
					strike_max = 0
					
					if BICs[len(BICs)-1] >= BICs[BIC_min]:
					#
						if len(BICs) > 2 and BICs[len(BICs)-2] <= BICs[len(BICs)-3]:
						#
							print("Carry min")
							print("Current BIC:")
							print(BIC_hmm)
							print("Current min:")
							print(BICs[BIC_min])
						#
						else:
						#
							print("Strike min")
							print("Current BIC:")
							print(BIC_hmm)
							print("Current min:")
							print(BICs[BIC_min])
							
							strike_min += 1
						#
					#
					else:
					#
						print("Reset min")
						print("Current BIC:")
						print(BIC_hmm)
						
						strike_min = 0
						
						BIC_min = len(BICs)-1
						
						best_hmm = deepcopy(hmm)
					#
				#
			#
			
			if n_clusters >= 10 and (strike_max >= 3 or strike_min >= 3):
			#
				print("Reached expected true minimum\n")
				
				break
			#
			
			print("\n")
		#
		
		if all_hmms[0] is None:
		#
			all_hmms[0] = deepcopy(best_hmm)
			all_BICs[0] = BICs
			all_BIC_mins[0] = BIC_min
			
			min_dim = n_dims
		#
		else:
		#
			all_hmms.append(deepcopy(best_hmm))
			all_BICs.append(BICs)
			all_BIC_mins.append(BIC_min)
			
			print("Current slowest timescale :")
			print(all_hmms[len(all_hmms)-1].timescales_[0])
			print("Previous slowest timescale :")
			print(all_hmms[len(all_hmms)-2].timescales_[0])
			
			if all_hmms[len(all_hmms)-1].timescales_[0] < all_hmms[len(all_hmms)-2].timescales_[0]:
			#
				print("Best HMM found with {:d} dimensions.".format(n_dims-1))
				
				min_dim = n_dims-1
				
				pos = len(all_hmms) - 2
				
				break
			#
			elif all_hmms[len(all_hmms)-1].timescales_[0] > all_hmms[len(all_hmms)-2].timescales_[0] and n_dims == n_comp:
			#
				print("Best HMM found with {:d} dimensions.".format(n_dims))
				
				min_dim = n_dims
				
				pos = len(all_hmms) - 1
				
				break
			#
		#
	#
	
	all_dims = np.arange(equil_dists, min_dim+2)
	
	if min_dim == 0:
	#
		raise Exception("Texte.")
	#
	
	return all_hmms, all_BICs, all_BIC_mins, min_dim, pos, all_dims
#
