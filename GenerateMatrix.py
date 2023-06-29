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
from sklearn.preprocessing import (
    MinMaxScaler,
    Normalizer,
    PowerTransformer,
    StandardScaler,
)
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from itertools import islice, product
from bisect import bisect, bisect_left, bisect_right
from collections import defaultdict

from functions import *

# use floor instead of comparison (smallest larger etc.)

def GenerateMatrix(
        chrom,
        CLAD,
        barcodes,
        df1,
        windows,
        peaks,
        peaks_file,
        hetsnp,
        snp_pos,
        window_size=1000000,
):
    windows_chr = windows.loc[windows["chrom"] == chrom][["chromStart", "chromEnd"]]
    windows_chr = windows_chr.astype('int')
    index = pd.MultiIndex.from_product(
        [df1["barcode"], windows_chr["chromStart"]],
        names=["barcode", "chromStart"],
    )
    result = pd.DataFrame(
        np.zeros(
            (
                len(df1["barcode"])
                * len(windows_chr["chromStart"]),
                3,
            )
        ),
        columns=["readcount", "H1", "H2"],
        index=index,
    )
    num_nonpeak = 0
    num_peak = 0
    waiting_reads = {}
    peak_seqnames = []
    peaks_chr = peaks.loc[peaks["chrom"] == chrom][["chromStart", "chromEnd"]]
    peaks_chr = peaks_chr.astype('int')
    lastwindow = 0

    for read in CLAD.fetch(chrom):  # parallel each chromosome
        if (
                not (read.is_unmapped) or not (read.mate_is_unmapped)
        ) and read.is_proper_pair:
            seqname = read.query_name
            if read.has_tag("CB"):
                barcode = read.get_tag("CB")
                if barcode not in barcodes:
                    continue

                startpos = read.reference_start
                endpos = read.reference_end
                chromosome = read.reference_name
                # Get belonging window
                startwindow = LargestSmaller(windows_chr.iloc[:, 0].tolist(), startpos)
                #startwindow = np.floor(startpos/window_size)
                lastpeak = LastPeak(peaks_chr, startpos)
                nextpeak = NextPeak(peaks_chr, startpos)
                closestsnp = SmallestLarger(snp_pos[chromosome], startpos)# Smallest SNP pos larger than startpos

                ######################################################################
                ############################ Window Ending ###########################
                ######################################################################

                # If proceeding to the next window
                if (
                        startwindow != lastwindow
                        and startwindow == lastwindow + window_size

                ):
                    print("Window_Ending", startwindow, lastwindow)
                    crosswindow_reads = {}
                    # Calculate the number of reads in peaks to resample
                    window = [chromosome, lastwindow, startwindow]
                    window = " ".join(str(x) for x in window)
                    chosenwindow = pybedtools.BedTool(window, from_string=True)
                    peak_frac = float(chosenwindow.coverage(peaks_file)[0][6])

                    # Calculating parameter for bernoulli distribution to randomly select reads in peak
                    try:
                        p = (
                                    (num_nonpeak / (1 - peak_frac) - num_nonpeak) / num_peak
                            ) / 2
                    except:
                        p = 0
                        print(num_nonpeak, peak_frac, num_peak, startwindow, lastwindow)

                    # Select reads in peaks
                    selected = np.array(peak_seqnames)[
                        np.array(bernoulli.rvs(p, size=len(peak_seqnames)), dtype=bool)
                    ]
                    selected = selected.tolist()
                    count = -1

                    # Loop through a copy and remove the pairs in waiting_reads
                    for current_seqname, current_read in waiting_reads.items():
                        count += 1

                        if len(current_read) == 1:  # Possibly pair in the next window
                            crosswindow_reads[current_seqname] = current_read
                            continue
                        else:
                            if current_seqname in selected:
                                for i in current_read:
                                    current_barcode = i[0]
                                    current_startwindow = i[3]
                                    
                                    # Readcount
                                    result.loc[
                                        (current_barcode, current_startwindow),
                                        "readcount",
                                    ] += 1

                                # Haplotype
                                haplotype_pair = ExtractNth(current_read, 5)
                                first_haplotype = current_read[0][5]
                                second_haplotype = current_read[1][5]
                                first_startwindow = current_read[0][3]
                                second_startwindow = current_read[1][3]
                                current_barcode = current_read[0][0]

                                if second_haplotype not in [
                                    -1,
                                    0,
                                ]:  # Pair has haplotype
                                    if (
                                            first_haplotype == haplotype_pair[1]
                                    ):  # If same haplotype and none false
                                        result.loc[
                                            (current_barcode, first_startwindow),
                                            "H" + str(first_haplotype),
                                        ] += 1
                                    elif first_haplotype in [-1, 0]:  # If current false
                                        result.loc[
                                            (current_barcode, second_startwindow),
                                            "H" + str(second_haplotype),
                                        ] += 1
                                    elif (
                                            second_haplotype != first_haplotype
                                            and first_haplotype in [1, 2]
                                    ):
                                        pass
                                else:  # Pair doesn't have haplotype
                                    if first_haplotype not in [-1, 0]:
                                        result.loc[
                                            (current_barcode, first_startwindow),
                                            "H" + str(first_haplotype),
                                        ] += 1
                                    else:
                                        pass

                    num_nonpeak = 0
                    num_peak = 0
                    waiting_reads = crosswindow_reads

                ######################################################################
                ############################### NON-PEAK #############################
                ######################################################################
                # Not in peak, excluding: start in previous peak, end in next peak, and covering whole peak(unlikely)
                if (
                        not startpos in range(lastpeak[0], lastpeak[1] + 1)
                        and not endpos in range(nextpeak[0], nextpeak[1] + 1)
                        and not (startpos <= nextpeak[0] and endpos >= nextpeak[1])
                ):
                    num_nonpeak += 1
                    # Store seqname of read to match non-peak and peak pairs
                    # waiting_reads: 'seqname':[[barcode, startpos, endpos, startwindow, peak (0/1), (haplotype)]]

                    ################### NON-Informative Read ####################
                    if (
                            closestsnp > endpos
                    ):  # When the read is not overlapping with a SNP but its pair might have haplotype
                        if seqname not in waiting_reads:  # Seqname has not been recorde
                            waiting_reads[seqname] = [
                                [barcode, startpos, endpos, startwindow, 0, -1]
                            ]
                        else:  # Seqname has been recorded
                            visited_read = waiting_reads[seqname][0]
                            visited_haplotype = visited_read[5]

                            # Record and skip if pair is in peak
                            if visited_read[4] == 1:
                                waiting_reads[seqname].append(
                                    [barcode, startpos, endpos, startwindow, 0, -1]
                                )
                            # If pair is not in peak
                            else:
                                # Haplotype
                                if visited_haplotype not in [
                                    0,
                                    -1,
                                ]:  # If pair has haplotype
                                    result.loc[
                                        (barcode, startwindow),
                                        "H" + str(visited_haplotype),
                                    ] += 1
                                else:  # If pair doesn't have haplotype
                                    pass

                                # Read count
                                result.loc[(barcode, startwindow), "readcount"] += 1
                                result.loc[(barcode, visited_read[3]), "readcount"] += 1

                                waiting_reads.pop(seqname)

                    ################# Informative Reads #################
                    else:  # When this read overlaps with a SNP
                        cigarList, seq = GenerateCIGARList(read)
                        haplotype = AssignHaplotype(read, hetsnp, cigarList, seq)

                        # If its mate have already been visited and recorded (need to add back those don't have pairs)
                        if (
                                seqname not in waiting_reads
                        ):  # Seqname has not been recorded
                            waiting_reads[seqname] = [
                                [
                                    barcode,
                                    startpos,
                                    endpos,
                                    startwindow,
                                    0,
                                    haplotype,
                                ]
                            ]
                        else:  # Seqname has been recorded
                            visited_read = waiting_reads[seqname][0]
                            visited_haplotype = visited_read[5]

                            # Record and skip if pair is in peak
                            if visited_read[4] == 1:
                                waiting_reads[seqname].append(
                                    [
                                        barcode,
                                        startpos,
                                        endpos,
                                        startwindow,
                                        0,
                                        haplotype,
                                    ]
                                )
                            # If pair is not in peak
                            else:
                                # Haplotype
                                if visited_haplotype not in [
                                    -1,
                                    0,
                                ]:  # Pair has haplotype
                                    if (
                                            visited_haplotype == haplotype
                                    ):  # If same haplotype and none false
                                        result.loc[
                                            (barcode, startwindow), "H" + str(haplotype)
                                        ] += 1
                                    elif haplotype in [-1, 0]:  # If current false
                                        result.loc[
                                            (barcode, visited_read[3]),
                                            "H" + str(visited_haplotype),
                                        ] += 1
                                    elif (
                                            visited_haplotype != haplotype
                                            and haplotype in [1, 2]
                                    ):
                                        pass
                                else:  # Pair doesn't have haplotype
                                    if haplotype not in [-1, 0]:
                                        result.loc[
                                            (barcode, startwindow), "H" + str(haplotype)
                                        ] += 1
                                    else:
                                        pass

                                # Read count
                                result.loc[(barcode, startwindow), "readcount"] += 1
                                result.loc[(barcode, visited_read[3]), "readcount"] += 1

                                waiting_reads.pop(seqname)

                ######################################################################
                ################################# PEAK ###############################
                ######################################################################
                else:
                    num_peak += 1
                    peak_seqnames.append(seqname)
                    if (
                            closestsnp > endpos
                    ):  # When the read is not overlapping with a SNP but its pair might have haplotype
                        if seqname in waiting_reads:
                            waiting_reads[seqname].append(
                                [barcode, startpos, endpos, startwindow, 1, -1]
                            )
                        else:
                            waiting_reads[seqname] = [
                                [barcode, startpos, endpos, startwindow, 1, -1]
                            ]
                    else:  # When this read overlaps with a SNP
                        cigarList, seq = GenerateCIGARList(read)
                        haplotype = AssignHaplotype(read, hetsnp, cigarList, seq)

                        if seqname in waiting_reads:
                            waiting_reads[seqname].append(
                                [
                                    barcode,
                                    startpos,
                                    endpos,
                                    startwindow,
                                    1,
                                    haplotype,
                                ]
                            )
                        else:
                            waiting_reads[seqname] = [
                                [
                                    barcode,
                                    startpos,
                                    endpos,
                                    startwindow,
                                    1,
                                    haplotype,
                                ]
                            ]

                lastwindow = startwindow
                
    startwindow = windows_chr.iloc[-1,-1]
    print("Window_Ending", startwindow, lastwindow)
    # Calculate the number of reads in peaks to resample
    window = [chromosome, lastwindow, startwindow]
    window = " ".join(str(x) for x in window)
    chosenwindow = pybedtools.BedTool(window, from_string=True)
    peak_frac = float(chosenwindow.coverage(peaks_file)[0][6])

    # Calculating parameter for bernoulli distribution to randomly select reads in peak
    try:
        p = (
                    (num_nonpeak / (1 - peak_frac) - num_nonpeak) / num_peak
            ) / 2
    except:
        p = 0
        print(num_nonpeak, peak_frac, num_peak, startwindow, lastwindow)

    # Select reads in peaks
    selected = np.array(peak_seqnames)[
        np.array(bernoulli.rvs(p, size=len(peak_seqnames)), dtype=bool)
    ]
    selected = selected.tolist()
    count = -1

    # Loop through a copy and remove the pairs in waiting_reads
    for current_seqname, current_read in waiting_reads.items():
        count += 1

        if len(current_read) == 1:  # Possibly pair in the next window
            continue
        else:
            if current_seqname in selected:
                for i in current_read:
                    current_barcode = i[0]
                    current_startwindow = i[3]

                    # Readcount
                    result.loc[
                        (current_barcode, current_startwindow),
                        "readcount",
                    ] += 1

                # Haplotype
                haplotype_pair = ExtractNth(current_read, 5)
                first_haplotype = current_read[0][5]
                second_haplotype = current_read[1][5]
                first_startwindow = current_read[0][3]
                second_startwindow = current_read[1][3]
                current_barcode = current_read[0][0]

                if second_haplotype not in [
                    -1,
                    0,
                ]:  # Pair has haplotype
                    if (
                            first_haplotype == haplotype_pair[1]
                    ):  # If same haplotype and none false
                        result.loc[
                            (current_barcode, first_startwindow),
                            "H" + str(first_haplotype),
                        ] += 1
                    elif first_haplotype in [-1, 0]:  # If current false
                        result.loc[
                            (current_barcode, second_startwindow),
                            "H" + str(second_haplotype),
                        ] += 1
                    elif (
                            second_haplotype != first_haplotype
                            and first_haplotype in [1, 2]
                    ):
                        pass
                else:  # Pair doesn't have haplotype
                    if first_haplotype not in [-1, 0]:
                        result.loc[
                            (current_barcode, first_startwindow),
                            "H" + str(first_haplotype),
                        ] += 1
                    else:
                        pass

    return result
