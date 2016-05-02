#!usr/bin/python

import os
import re
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import numpy as np
from argparse import ArgumentParser
from copy import deepcopy

parser = ArgumentParser()

parser.add_argument("-i1", "--input1", nargs='*') # .. # Files in the input1 argument will be printed in separate plots
parser.add_argument("-i2", "--input2", nargs='*') # .. # Files in the input2 argument will be printed in every plot
parser.add_argument("-o", "--output")

args = vars(parser.parse_args())

input1 = args["input1"]
input2 = args["input2"]
plt_out = args["output"]

if not input1 is None:
#
	for file1 in input1 :
	#
		if not os.path.isfile(file1):
		#
			raise Exception("File {:s} does not exist.".format(file1))
		#
	#
#

if not input2 is None:
#
	for file2 in input2 :
	#
		if not os.path.isfile(file2):
		#
			raise Exception("File {:s} does not exist.".format(file2))
		#
	#
#

formats_2 = []
xnames_2 = []
x_2 = []
ynames_2 = []
y_2 = []

for file2 in input2:
#
	with open(file2, 'r') as rf:
	#
		x = []
		ynames = []
		y = []
		y2 = []
		
		formats = []
		format2 = "line"
		
		datatype = "x"
		
		for line in rf:
		#
			if re.match("\d", line):
			#
				if datatype == "y":
				#
					y2.append(float(line))
				#
				elif datatype == "x":
				#
					x.append(float(line))
				#
			#
			elif re.match("format:", line):
			#
				format2 = re.match("format:(\w+)", line).group(1)
			#
			elif re.match("[^X0-9]", line):
			#
				if x:
				#
					x_c = deepcopy(x)
					
					x_2.append(x_c)
					
					x = []
				#
				
				if y2:
				#
					y2_c = deepcopy(y2)
					
					y.append(y2_c)
					
					y2 = []
				#
				
				format_c = deepcopy(format2)
				
				formats.append(format_c)
				
				format2 = "line"
				
				dataname = deepcopy(line[:-1])
				
				ynames.append(dataname)
				
				datatype = "y"
			#
			elif re.match("X:", line):
			#
				xnames_2.append(re.match("X:([\w\_]+)", line).group(1))
				
				if y:
				#
					y_c = deepcopy(y)
					
					y_2.append(y_c)
					
					y = []
					
					ynames_c = deepcopy(ynames)
					
					ynames_2.append(ynames_c)
					
					ynames = []
					
					formats_c = deepcopy(formats)
					
					formats_2.append(formats_c)
					
					formats = []
				#
				
				datatype = "x"
			#
		#
		
		y.append(y2)
		
		print(y)
		
		y_c = deepcopy(y)
		
		y_2.append(y_c)
		
		print(ynames)
		
		formats_c = deepcopy(formats)
					
		formats_2.append(formats_c)
		
		ynames_c = deepcopy(ynames)
		
		ynames_2.append(ynames_c)
	#
#


colors = [cm.jet(float(i)/5.0) for i in range(6)]

barstat = 0

pointstat = 0

for i in range(len(x_2)):
#
	for j in range(len(y_2[i])):
	#
		ymax = 0
		
		ymin = math.fabs(y_2[i][j][0])
		
		for k in range(len(y_2[i][j])):
		#
			if math.fabs(y_2[i][j][k]) > ymax:
			#
				ymax = math.fabs(y_2[i][j][k])
			#
			
			if math.fabs(y_2[i][j][k]) < ymin:
			#
				ymin = math.fabs(y_2[i][j][k])
			# 
		#
		
		for k in range(len(y_2[i][j])):
		#
			y_2[i][j][k] = (y_2[i][j][k] - ymin)/(ymax-ymin)
		# 
		
		if formats_2[i][j] == "line":
		#
			plt.plot(x_2[i], y_2[i][j], label=ynames_2[i][j])
		#
		elif formats_2[i][j] == "bar":
		#
			x_temp = deepcopy(x_2[i])
			
			for k in range(len(x_temp)):
			#
				x_temp[k] -= 0.5
			#
			
			print(x_temp)
			
			print(y_2[i][j])
			
			print(ynames_2[i][j])
			
			plt.bar(x_temp, y_2[i][j], label=ynames_2[i][j], color = colors[barstat])
			
			barstat += 1
		#
	#
#

plt.title("Relationship between tICA scans and Kex from CPMG\nfor the 7 relevant eigenvectors.")

plt.xlabel("Residue")

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.show()
