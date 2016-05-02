#!usr/bin/python

from argparse import ArgumentParser
import os
from checkpoint.io_results import load_feat, load_tica, load_hmm
import numpy as np
import math

parser = ArgumentParser()

parser.add_argument("-f1", "--feat1")
parser.add_argument("-f2", "--feat2", nargs='*')
parser.add_argument("-t1", "--tica1")
parser.add_argument("-t2", "--tica2", nargs='*')
#parser.add_argument("-h1", "--hmm1")
#parser.add_argument("-h2", "--hmm2")
#parser.add_argument("-o", "--output")

args = vars(parser.parse_args())

del parser

feat1 = args["feat1"]
feats2 = args["feat2"]
tica1 = args["tica1"]
ticas2 = args["tica2"]
#hmm1 = args["hmm1"]
#hmm2 = args["hmm2"]

if not os.path.isfile(feat1):
#
	raise Exception("tICA {:s} does not exist.".format(feat1))
#

for feat2 in feats2 :
#
	if not os.path.isfile(feat2):
	#
		raise Exception("tICA {:s} does not exist.".format(feat2))
	#
#

if not os.path.isfile(tica1):
#
	raise Exception("tICA {:s} does not exist.".format(tica1))
#

for tica2 in ticas2 :
#
	if not os.path.isfile(tica2):
	#
		raise Exception("tICA {:s} does not exist.".format(tica2))
	#
#

#if not os.path.isfile(hmm1):
#
	#raise Exception("tICA {:s} does not exist.".format(hmm1))
#

#if not os.path.isfile(hmm2):
#
	#raise Exception("tICA {:s} does not exist.".format(hmm2))
#

ft1 = load_feat(feat1)
td1 = load_tica(tica1)

pre_evecs1 = np.asarray(td1['eigenvectors']).T.tolist()

evecs1_full = pre_evecs1[:td1['n_comp']]
evals1_full = td1['eigenvalues'][:td1['n_comp']]

for k in range(len(ticas2)):
#
	print("Loading {:s} with {:s} ...".format(ticas2[k], feats2[k]))
	
	ft2 = load_feat(feats2[k])
	td2 = load_tica(ticas2[k])
	
	pre_evecs2 = np.asarray(td2['eigenvectors']).T.tolist()

	evecs2_full = pre_evecs2[:td2['n_comp']]
	evals2_full = td2['eigenvalues'][:td2['n_comp']]

	matches = dict()

	for i in range(len(ft1)):
	#
		for j in range(len(ft2)):
		#
			if ft1[i] == ft2[j]:
			#
				matches[i] = j
			
				break
			#
		#
	#

	evecs1 = [[0.0 for j in range(len(matches))] for i in range(len(evecs1_full))]
	evecs2 = [[0.0 for j in range(len(matches))] for i in range(len(evecs2_full))]

	for i in range(len(evecs1)):
	#
		j = -1
	
		for match in matches:
		#
			j += 1
		
			evecs1[i][j] = evecs1_full[i][match]
		#
	#

	for i in range(len(evecs2)):
	#
		j = -1
	
		for match in matches:
		#
			j += 1
		
			evecs2[i][j] = evecs2_full[i][matches[match]]
		#
	#
	
	print("Comparing {:s} with {:s}".format(tica1, ticas2[k]))

	for i in range(len(evecs1)):
	#
		for j in range(len(evecs2)):
		#
			score = np.vdot(evecs1[i], evecs2[j])/math.sqrt(np.vdot(evecs1[i], evecs1[i]) * np.vdot(evecs2[j], evecs2[j]))
		
			if (score > math.sqrt(2.0)/2.0) or (score < -math.sqrt(2.0)/2.0):
			#
				print("{:d}\t{:d}\t{:f}".format(i, j, score))
			#
		#
	#
#
