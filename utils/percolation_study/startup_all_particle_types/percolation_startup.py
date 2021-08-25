import os 
import argparse


#TODO:? command line tool below plus function tool 

parser = argparse.ArgumentParser()
parser.add_argument('-ptype', type=str, choices=['dma-as2', 'feq-as2', 'dma-as1', 'dmo-s1', 'dmo-s2', 'dmo-as1'])
parser.add_argument('-delta', type=float)
parser.add_argument('-radius', type=float)
parser.add_argument('-source-dest',type=str)
parser.add_argument('-new-dir-name', type=str)





args = parser.parse_args()
