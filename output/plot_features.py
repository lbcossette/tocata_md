import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
import math
from matplotlib.patches import Ellipse


########################################################################
#																	   #
#	 			Plot calculated slopes for each component			   #
#																	   #
########################################################################


def plot_slopes(traj_tica, max_i, savepath):
#
	slopes = None
	
	for i in range(2, max_i - 2):
	#
		local_slope = (traj_tica.timescales_[i+2] - traj_tica.timescales_[i-2])/8.0 + (traj_tica.timescales_[i+2] + traj_tica.timescales_[i+1] - traj_tica.timescales_[i-1] - traj_tica.timescales_[i-2])/12.0
		
		if slopes is None:
		#
			slopes = np.array([local_slope])
		#
		else:
		#
			slopes = np.append(slopes, [local_slope])
		#
	#
	
	plt.plot(np.arange(2, max_i - 2, 1), slopes)
	
	plt.savefig(savepath)
	plt.close()
#


########################################################################
#																	   #
#	 		   Plot calculated timescales for each component		   #
#																	   #
########################################################################


def plot_timescales(traj_tica, max_i, savepath):
#
	timescales = None
	
	for i in range(0, max_i):
	#
		if timescales is None:
		#
			timescales = np.array([traj_tica.timescales_[i]])
		#
		else:
		#
			timescales = np.append(timescales, [traj_tica.timescales_[i]])
		#
	#
	
	plt.plot(np.arange(0, max_i, 1), timescales)
	
	plt.savefig(savepath)
	plt.close()
#


########################################################################
#																	   #
#	 	   Get a manageable number of points through splitting		   #
#																	   #
########################################################################


def manageable(length):
#
	return int(length/90000 + 1)
#


########################################################################
#																	   #
#	 	   		  Get the estimated range of a dataset				   #
#																	   #
########################################################################


def get_range(means, sigmas):
#
	extrema = [[m+6.0*v for m,v in zip(means[0], sigmas[0])], [m-6.0*v for m,v in zip(means[0], sigmas[0])]]
	
	for i in range(len(means)):
	#
		extrema.append([m+6.0*v for m,v in zip(means[i], sigmas[i])])
		extrema.append([m-6.0*v for m,v in zip(means[i], sigmas[i])])
	#
	
	x_range = [extrema[0][0], extrema[0][0]]
	y_range = [extrema[0][1], extrema[0][1]]
	
	for i in range(len(extrema)):
	#
		if extrema[i][0] > x_range[1]:
		#
			x_range[1] = extrema[i][0]
		#
		elif extrema[i][0] < x_range[0]:
		#
			x_range[0] = extrema[i][0]
		#
		
		if extrema[i][1] > y_range[1]:
		#
			y_range[1] = extrema[i][1]
		#
		elif extrema[i][1] < y_range[0]:
		#
			y_range[0] = extrema[i][1]
		#
	#
	
	return x_range, y_range
#


########################################################################
#																	   #
#	 		  Plot point density for each pair of components		   #
#																	   #
########################################################################


def plot_density(ic_slice, skip, savepath):
#
	print("skip :")
	print(skip)
	
	all_trajs = ic_slice[0]
	
	print(np.shape(all_trajs))
	
	for i in range(1, len(ic_slice)):
	#
		all_trajs = np.append(all_trajs, ic_slice[i], axis=0)
		
		#print(np.shape(all_trajs))
	#
	
	all_trajs = np.split(all_trajs, 2, axis=1)
	
	print(np.shape(all_trajs))
	
	points_x = np.array([all_trajs[0][0][0]])
	
	for i in range(1, len(all_trajs[0])):
	#
		if i % skip == 0:
		#
			points_x = np.append(points_x, [all_trajs[0][i][0]])
		#
	#
	
	print("points_x :")
	print(np.shape(points_x))
	
	points_y = np.array([all_trajs[1][0][0]])
	
	for i in range(1, len(all_trajs[0])):
	#
		if i % skip == 0:
		#
			points_y = np.append(points_y, [all_trajs[1][i][0]])
		#
	#
	
	print("points_y :")
	print(np.shape(points_y))
	
	del all_trajs
	
	print("Stacking points...")
	
	points_xy = np.vstack([points_x,points_y])
	
	print("Calculating point density...")
	
	points_z = gaussian_kde(points_xy)(points_xy)
	
	print("Sorting points...")
	
	idx = points_z.argsort()
	
	x, y, z = points_x[idx], points_y[idx], points_z[idx]
	
	print("Rendering...")
	
	fig, ax = plt.subplots()
	
	ax.scatter(x, y, c=z, s=5, edgecolor='')
	
	plt.savefig(savepath)
	plt.close()
#


########################################################################
#																	   #
#	 		  				Plot HMM clusters						   #
#																	   #
########################################################################


def plot_hmm(ic_slice, means, sigmas, skip, savepath):
#
	print("Generating ellipses...")
	
	ells = [Ellipse(xy = [means[i][0], means[i][1]], width=1.96*math.sqrt(sigmas[i][0]), height=1.96*math.sqrt(sigmas[i][1]), angle=0) for i in range(len(means))]
	
	fig, ax = plt.subplots()
	
	for i in range(len(ells)):
	#
		ax.add_artist(ells[i])
		ells[i].set_clip_box(ax.bbox)
		ells[i].set_alpha(0.3)
		ells[i].set_facecolor('black')
	#
	
	print("Generating scatterplot...")
	
	print("skip :")
	print(skip)
	
	all_trajs = ic_slice[0]
	
	print(np.shape(all_trajs))
	
	for i in range(1, len(ic_slice)):
	#
		all_trajs = np.append(all_trajs, ic_slice[i], axis=0)
		
		#print(np.shape(all_trajs))
	#
	
	all_trajs = np.split(all_trajs, 2, axis=1)
	
	print(np.shape(all_trajs))
	
	points_x = np.array([all_trajs[0][0][0]])
	
	for i in range(1, len(all_trajs[0])):
	#
		if i % skip == 0:
		#
			points_x = np.append(points_x, [all_trajs[0][i][0]])
		#
	#
	
	print("points_x :")
	print(np.shape(points_x))
	
	points_y = np.array([all_trajs[1][0][0]])
	
	for i in range(1, len(all_trajs[0])):
	#
		if i % skip == 0:
		#
			points_y = np.append(points_y, [all_trajs[1][i][0]])
		#
	#
	
	print("points_y :")
	print(np.shape(points_y))
	
	del all_trajs
	
	print("Stacking points...")
	
	points_xy = np.vstack([points_x,points_y])
	
	print("Calculating point density...")
	
	points_z = gaussian_kde(points_xy)(points_xy)
	
	print("Sorting points...")
	
	idx = points_z.argsort()
	
	x, y, z = points_x[idx], points_y[idx], points_z[idx]
	
	ax.scatter(x, y, c=z, s=5, edgecolor='')
	
	print("Rendering...")
	
	plt.savefig(savepath)
	plt.close()
#


########################################################################
#																	   #
#	 		  	Evaluate the density function of a HMM				   #
#																	   #
########################################################################


def density(x, y, means, sigma, weights):
#
	value = 0.0
	
	for i in range(len(means)):
	#
		arg_exp = 0.5*((x - means[i][0])**2 / sigma[i][0] + (y - means[i][1])**2 / sigma[i][1])
		
		divisor = 2 * np.pi * np.sqrt(sigma[i][0] * sigma[i][1])
		
		value += weights[i] * np.exp(-1.0 * arg_exp) / divisor
	#
	
	return value
#


########################################################################
#																	   #
#	 		  		   Plot HMM density estimate					   #
#																	   #
########################################################################


def plot_estimate(means, sigma, weights, savepath):
#
	print("Generating density estimate...")
	
	image = np.empty([1000,1000])
	
	x_range, y_range = get_range(means, sigma)
	
	for i in range(1000):
	#
		for j in range(1000):
		#
			x1 = x_range[0] + i*(x_range[1] - x_range[0])/1000
			x2 = x_range[0] + (i+1)*(x_range[1] - x_range[0])/1000
			y1 = y_range[0] + j*(y_range[1] - y_range[0])/1000
			y2 = y_range[0] + (j+1)*(y_range[1] - y_range[0])/1000
			
			dist_avg = density(x1, y1, means, sigma, weights) + density(x2, y1, means, sigma, weights) + density(x1, y2, means, sigma, weights) + density(x2, y2, means, sigma, weights)
			
			dist_avg /= 4.0
			
			image[i][j] = dist_avg
		#
	#
	
	plt.imsave(savepath, image, cmap="jet")
	plt.close()
#


########################################################################
#																	   #
#	 		  		     Plot population counts 					   #
#																	   #
########################################################################


def plot_popl_1d(ic_slice, savepath, dimension):
#
	all_trajs = ic_slice[0]
	
	for i in range(1, len(ic_slice)):
	#
		all_trajs = np.append(all_trajs, ic_slice[i], axis=0)
	#
	
	x_range = [all_trajs[0][0],all_trajs[0][0]]
	
	for i in range(len(all_trajs)):
	#
		if x_range[0] > all_trajs[i][0]:
		#
			x_range[0] = all_trajs[i][0]
		#
		
		if x_range[1] < all_trajs[i][0]:
		#
			x_range[1] = all_trajs[i][0]
		#
	#
	
	values = np.zeros(dimension)
	pos_values = np.zeros(dimension)
	
	for i in range(len(all_trajs)):
	#
		x_pos = int((all_trajs[i][0] - x_range[0])/(x_range[1] - x_range[0]) * float(dimension))
		
		if x_pos == dimension:
		#
			x_pos -= 1
		#
		
		values[x_pos] += 1
	#
	
	for i in range(len(values)):
	#
		values[i] /= float(len(all_trajs))
		
		pos_values[i] = (float(i) + 0.5)*(x_range[1] - x_range[0])/float(dimension)
	#
	
	plt.plot(pos_values, values)
	
	plt.savefig(savepath)
	plt.close()
#


########################################################################
#																	   #
#	 		   Evaluate the density function of a 1D HMM			   #
#																	   #
########################################################################


def density_1d(x, means, sigma, weights):
#
	value = 0.0
	
	for i in range(len(means)):
	#
		arg_exp = 0.5*((x - means[i][0])**2 / sigma[i][0])
		
		divisor = 2 * np.pi * np.sqrt(sigma[i][0])
		
		value += weights[i] * np.exp(-1.0 * arg_exp) / divisor
	#
	
	return value
#


########################################################################
#																	   #
#	 		  		   Plot HMM 1D density estimate					   #
#																	   #
########################################################################


def plot_estimate_1d(means, sigma, weights, x_range, savepath, dimension):
#
	print("Generating density estimate...")
	
	values = np.zeros(dimension)
	pos_values = np.zeros(dimension)
	
	for i in range(dimension):
	#
		x1 = x_range[0] + i*(x_range[1] - x_range[0])/dimension
		x2 = x_range[0] + (i+1)*(x_range[1] - x_range[0])/dimension
		
		dist_avg = density_1d(x1, means, sigma, weights) + density_1d(x2, means, sigma, weights)
		
		dist_avg /= 2.0
		
		values[i] = dist_avg
		
		pos_values[i] = (float(i) + 0.5)*(x_range[1] - x_range[0])/float(dimension)
	#
	
	plt.plot(pos_values, values)
	
	plt.savefig(savepath)
	plt.close()
#


########################################################################
#																	   #
#	 		  		Plot BICs vs number of clusters					   #
#																	   #
########################################################################


def plot_BICs(all_BICs, all_dims, savepath):
#
	print("Plotting BICs...")
	
	for i in range(len(all_BICs)):
	#
		n_clusters = np.arrange(1,len(all_BICs[i])+1)
		
		plt.plot(n_clusters, all_BICs[i], label = "{:d} dimensions".format(all_dims[i]))
	#
	
	plt.savefig(savepath)
	plt.close()
#
