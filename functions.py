import numpy as np
import random
import pybedtools
import sys
import os
from IPython.utils import io
import pysam
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import seaborn as sns
import math
from scipy.stats import zscore, bernoulli
from sklearn.preprocessing import MinMaxScaler, Normalizer, PowerTransformer
from sklearn.preprocessing import Normalizer
from sklearn.preprocessing import PowerTransformer
from sklearn.preprocessing import StandardScaler
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from itertools import islice, product
from bisect import bisect, bisect_left, bisect_right
from collections import defaultdict


def InRange(mylist, pos):
    for idx, val in enumerate(mylist):
        if(pos in range(val[0],val[1]+1)):
            return True
    return False

def WhichRange(mylist, pos):
    diff = np.zeros(len(mylist))
    for idx, val in enumerate(mylist):
        diff[idx] = val[1] - val[0] + 1
    for idx, val in enumerate(mylist):
        if(pos in range(val[0],val[1]+1)):
            return [idx, sum(islice(diff, idx))]
    return -1

def LargestSmaller(myList, pos):
    i = bisect(myList, pos)
    if i:
        return myList[i-1]
    else:
        return -1

def SmallestLarger(myList, pos):
    i = bisect(myList, pos)
    if i >= len(myList):
        return np.inf
    elif i:
        return myList[i]
    else:
        return -1

def LastPeak(peaks, pos):
    start = peaks.iloc[:,0].to_numpy()
    idx = bisect(start, pos) - 1
    if idx != -1:
        return peaks.iloc[idx].values.flatten().tolist()
    else:
        return [0,1]

def NextPeak(peaks, pos):
    start = peaks.iloc[:,0].to_numpy()
    idx = bisect(start, pos)
    if idx < len(start):
        return peaks.iloc[idx].values.flatten().tolist()
    else:
        return [999999998,999999999]

def WindowPercent(cigarList, seq, window_size):
    length = len(seq)
    lastpos = cigarList[len(cigarList)-1][1]
    firstpos = cigarList[0][0]
    lastwindow = lastpos // window_size
    firstwindow = firstpos // window_size

    windows = []
    for i in range(firstwindow,lastwindow+1):
        window_start = i * window_size
        window_end = window_start + window_size - 1

        read_coverage = 0
        for pair in cigarList:
            if pair[1] <= window_start:
                continue
            if pair[0] >= window_end:
                break

            overlap_start = max(pair[0], window_start)
            overlap_end = min(pair[1], window_end)
            read_coverage += overlap_end - overlap_start + 1

        coverage_percentage = read_coverage / length
        windows.append((window_start, coverage_percentage))

    return windows

def GenerateCIGARList(read):
    cigarLine = read.cigar;
    cigarList = []
    cur_pos = read.reference_start + 1
    seq = read.seq
    if len(cigarLine) == 1 and cigarLine[0][0] == 0:
        pass
    else:
        for (cigarType,cigarLength) in cigarLine:
            try:
                if(cigarType == 0): # Match
                    next_pos = cur_pos + cigarLength - 1
                    cigarList.append([cur_pos, next_pos])
                    cur_pos = next_pos + 1
                elif(cigarType == 1): # Insertions # Should I remove it?
                    ins_pos = cur_pos - (read.reference_start + 1)
                    seq = seq[:ins_pos] + seq[ins_pos+cigarLength:]
                elif(cigarType == 2): # Deletion
                    cur_pos = cur_pos + cigarLength
                elif(cigarType == 3): # Skip
                    cur_pos = cur_pos + cigarLength
                elif(cigarType == 4): # Soft clipping
                    sc_pos = cur_pos - (read.reference_start + 1)
                    seq = seq[:sc_pos] + seq[sc_pos+cigarLength:]
                elif(cigarType == 5): # Hard clipping
                    continue
                elif(cigarType == 6): # Padding
                    continue
                else: # Catch and expect later?
                    print("Wrong CIGAR number")
                    sys.exit(1);
            except:
                print("Problem")
                sys.exit(1);
    return cigarList, seq

def AssignHaplotype(read, hetsnp, cigarList, seq):
    haplotype = 0
    for snp in hetsnp.fetch(read.reference_name, read.reference_start, read.reference_end):
        # If SNP is in gap of read (N or D)
        if len(cigarList) != 0 and InRange(cigarList, snp.pos) == False:
            continue
        # If only matches or in range of read
        else:
            gts = [s['GT'] for s in snp.samples.values()]
            if len(cigarList) == 0:
                startpos = read.reference_start + 1
                diffsum = 0
            else:
                # Start of other ranges
                idx, diffsum = WhichRange(cigarList, snp.pos)
                startpos = cigarList[idx][0]
            # Get the base in read at the SNP position
            base = seq[snp.pos - startpos + int(diffsum)]

            # If read have not been visited (first SNP overlapping this read)
            if haplotype == 0:
                if gts == [(0, 1)]:
                    if base == snp.alleles[0]:
                        haplotype = 1
                    elif base == snp.alleles[1]:
                        haplotype = 2
                    else: # No match
                        pass
                elif gts == [(1, 0)]:
                    if base == snp.alleles[1]:
                        haplotype = 1
                    elif base == snp.alleles[0]:
                        haplotype = 2
                    else: # No match
                        pass
            else: # If multiple SNPs overlapping with this read
                if gts == [(0, 1)]:
                    # If same haplotype then do nothing
                    if base == snp.alleles[0] and haplotype == 1:
                        pass
                    elif base == snp.alleles[1] and haplotype == 2:
                        pass
                    elif base != snp.alleles[0] or base != snp.alleles[1]:
                        pass
                    else: # If a read overlaps with two different haplotypes, ignore this read
                        haplotype = 0
                        break
                elif gts == [(1, 0)]:
                    if base == snp.alleles[1] and haplotype == 1:
                        pass
                    elif base == snp.alleles[0] and haplotype == 2:
                        pass
                    elif base != snp.alleles[0] or base != snp.alleles[1]: # if miss match, ignore this snp
                        pass
                    else: # If a read overlaps with two different haplotypes, ignore this read
                        haplotype = 0
                        break
    return haplotype

def ExtractNth(myList, i):
    return [item[i] for item in myList]