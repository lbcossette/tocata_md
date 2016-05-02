#!usr/bin/python

from analysis.inheritance import inheritance
from checkpoint.io_results import load_hmm
import os
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument("-i", "--input")

args = vars(parser.parse_args())

del parser

hmm_dict = load_hmm(args["input"])

tree = inheritance(hmm_dict['means'], hmm_dict['vars'], hmm_dict['popls'], hmm_dict['transmat'])

tree.lump()

tree.retrolump()
