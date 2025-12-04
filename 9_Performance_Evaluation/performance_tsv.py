from sys import argv
import math
import os

'''
performance_tsv.py

This script computes the confusion matrix, MCC, tpr, ppv and generates a .tsv file
as an output for further data analysis. It takes as input a .class file.
*The treshold values are set in the function run_performance().

Usage:
    python performance_tsv.py <set_#.class>
'''

# The first part is the same as the performance.py script

# Get the confusion matrix
def get_cm(filename, threshold, pe, pr=1):
    cm = [[0, 0], [0, 0]]
    with open(filename) as f:
        for line in f:
            v = line.rstrip().split()
            evalue = float(v[pe])
            r = int(v[pr])
            p = 1 if evalue <= threshold else 0
            cm[p][r] += 1
    return cm

# Compute accuracy
def get_q2(cm):
    n = sum([sum(row) for row in cm])
    return (cm[0][0] + cm[1][1]) / n if n > 0 else 0

# Compute Matthews Correlation Coefficient
def get_mcc(cm):
    tn, fn = cm[0]
    fp, tp = cm[1]
    numerator = tp * tn - fp * fn
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return numerator / denominator if denominator != 0 else 0

# True positive rate
def get_tpr(cm):
    tp, fn = cm[1][1], cm[1][0]
    return tp / (tp + fn) if (tp + fn) != 0 else 0

# Positive predictive value (precision)
def get_ppv(cm):
    tp, fp = cm[1][1], cm[0][1]
    return tp / (tp + fp) if (tp + fp) != 0 else 0

# Run the script using different tresholds. Iterates from 1e-(1 to 12).
def run_performance(filename, fullseq=True):
    pe = 2 if fullseq else 3  # Column index for E-value
    results = []
    for i in range(1, 13):  # 1e-1 to 1e-12
        threshold = float(f"1e-{i}")
        cm = get_cm(filename, threshold, pe)
        q2 = get_q2(cm)
        mcc = get_mcc(cm)
        tpr = get_tpr(cm)
        ppv = get_ppv(cm)
        results.append([f"1e-{i}", q2, mcc, tpr, ppv, fullseq])
    return results

def write_results_to_tsv(filename, results):
    basename = os.path.splitext(filename)[0]
    output_file = f"{basename}_performance.tsv"
    with open(output_file, 'w') as out:
        out.write("threshold\tq2\tmcc\ttpr\tppv\tfullseq\n")
        for row in results:
            out.write("\t".join(map(str, row)) + "\n")
    print(f"Results written to {output_file}")

if __name__ == '__main__':
    if len(argv) < 2:
        print("ERROR: This program takes a .class file as an argument.")
        print("Usage: python performance.py <class_file>")
        exit(1)

    class_file = argv[1]
    results_full = run_performance(class_file, fullseq=True)
    results_dom = run_performance(class_file, fullseq=False)
    write_results_to_tsv(class_file, results_full + results_dom)
