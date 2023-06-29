import numpy as np
from multiprocessing import Process, Manager
import concurrent.futures
import pickle
import sys
from GenerateMatrix import *

if __name__ == "__main__":
    folder = '/scratch/remills_root/remills0/wjhlang/cnv_caller/CLAD/'

    CLAD = os.path.join(folder, 'CLAD_rmdup_sorted.bam')
    CLAD = pysam.AlignmentFile(CLAD, mode='rb')

    peaks_path = os.path.join(folder, 'atac_peaks_chrom.bed')
    peaks_file = pybedtools.BedTool(peaks_path)

    windows = pd.read_csv(os.path.join(folder, 'windows_1m.bed'), sep='\t')
    windows = windows.astype({"chromStart": int, "chromEnd": int})


    CB = list()
    for read in CLAD.fetch():
        if read.has_tag('CB'):
            CB.append(read.get_tag('CB'))

    CB = np.array(CB)
    unique, counts = np.unique(CB, return_counts=True)
    df = pd.DataFrame({"barcode": unique, "counts": counts})
    df = df.sort_values(by='counts', ascending=False)
    df1 = df.head(6000)
    barcodes = list(df1["barcode"])

    # Peaks
    peaks = []
    with open(peaks_path) as f:
        for line in f:
            x = line.strip().split()
            if len(x) == 3:
                peaks.append(x)
    peaks = pd.DataFrame(peaks, columns=['chrom', 'chromStart', 'chromEnd'])
    peaks = peaks.astype({"chromStart": int, "chromEnd": int})

    # SNPs
    hetsnp = os.path.join(folder, 'CLAD_hetsnp.vcf.gz')
    hetsnp = pysam.VariantFile(hetsnp)

    snp_pos = defaultdict(list)
    for snp in hetsnp.fetch():
        snp_pos[snp.chrom].append(snp.pos)

    # Result
    result_all = {}

    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7','chr8','chr9', 'chr10',
              'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
              'chr20', 'chr21', 'chr22']

    # def run_matrix(chrom):
    #     return chrom, GenerateMatrix(chrom, CLAD, barcodes, df1, windows, peaks, peaks_file, hetsnp, snp_pos)
    #
    # result_all = defaultdict(list)
    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     futures = {executor.submit(run_matrix, chrom): chrom for chrom in chroms}
    #     for future in concurrent.futures.as_completed(futures):
    #         try:
    #             k, v = future.result()
    #         except Exception as e:
    #             print(f"{futures[future]} throws {e}")
    #         else:
    #             result_all[k] = v
    #
    # with open('result_all.pkl', 'wb') as f:
    #     pickle.dump(result_all, f)
    chrom = sys.argv[1]
    result = GenerateMatrix(chrom, CLAD, barcodes, df1, windows, peaks, peaks_file, hetsnp, snp_pos)
    result.to_csv('result' + chrom + '.csv')
    # with open('result_chr10.pkl', 'wb') as f:
    #  pickle.dump(result, f)

    # def run_matrix(chrom):
    #     return chrom, GenerateMatrix(chrom, CLAD, barcodes, df1, windows, peaks, peaks_file, hetsnp, snp_pos)
    #
    # result_all = defaultdict(list)
    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     futures = {executor.submit(run_matrix, chrom): chrom for chrom in chroms}
    #     for future in concurrent.futures.as_completed(futures):
    #         try:
    #             chrom, result = future.result()
    #         except Exception as e:
    #             print(f"{futures[future]} throws {e}")
    #         else:
    #             result.to_csv(chrom + '.csv')

