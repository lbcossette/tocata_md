import numpy as np
import math





def compare_ticas(evecs1, tscales1, features1, evecs2, tscales2, features2):
#
	matches = dict()

	for i in range(len(features1)):
	#
		for j in range(len(features2)):
		#
			if features1[i] == features2[j]:
			#
				matches[i] = j
			
				break
			#
		#
	#
	
	evecs_trim1 = [[0.0 for j in range(len(matches))] for i in range(len(evecs1))]
	evecs_trim2 = [[0.0 for j in range(len(matches))] for i in range(len(evecs2))]

	for i in range(len(evecs1)):
	#
		j = -1
	
		for match in matches:
		#
			j += 1
		
			evecs_trim1[i][j] = evecs1[i][match]
		#
	#

	for i in range(len(evecs2)):
	#
		j = -1
	
		for match in matches:
		#
			j += 1
		
			evecs_trim2[i][j] = evecs2[i][matches[match]]
		#
	#
	
	coss = [0.0 for i in range(len(evecs_trim1))]
	rtimes = [0.0 for i in range(len(evecs_trim1))]
	
	for i in range(len(evecs_trim1)):
	#
		for j in range(len(evecs_trim2)):
		#
			score = np.vdot(evecs_trim1[i], evecs_trim2[j])/math.sqrt(np.vdot(evecs_trim1[i], evecs_trim1[i]) * np.vdot(evecs_trim2[j], evecs_trim2[j]))
			
			if math.fabs(score) > math.fabs(coss[i]):
			#
				coss[i] = score
				
				rtimes[i] = tscales2[j]/tscales1[i]
			#
		#
	#
	
	return coss, rtimes
#
