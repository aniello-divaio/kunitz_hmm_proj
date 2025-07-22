from sys import argv
import math

'''
performance.py

This script computes the confusion matrix, MCC, tpr, ppv. It takes as input a
.class file and a treshold value as arguments.

A third parameter can be set if performance evaluation of full-sequence, single domain or both is desired:
1 -> full-sequence
2 -> single domain
0 or undefined -> both

Usage:
    python performance.py <set_#.class> <1e-$i>
'''

# Get the confusion matrix
def get_cm(filename, threshold, pe, pr=1):
    cm = [[0, 0], [0, 0]]
    f = open(filename)
    for line in f:
        v = line.rstrip().split()
        # v[pe] is the E-value, v[pr] is the true class (1 = Kunitz, 0 = non-Kunitz)
        evalue = float(v[pe])
        r = int(v[pr])
        # Predict class based on threshold: below = Kunitz (1), above = non-Kunitz (0)
        if evalue <= threshold:
            p = 1
        else:
            p = 0
        # Fill confusion matrix (prediction, real class)
        cm[p][r] = cm[p][r] + 1
    return cm

# Compute accuracy
def get_q2(cm):
    n = float(cm[0][0] + cm[0][1] + cm[1][0] + cm[1][1])
    return (cm[0][0] + cm[1][1])/n

# Compute Matthews Correlation Coefficient
def get_mcc(cm):
    d = math.sqrt(
        (cm[0][0] + cm[1][0]) *
        (cm[0][0] + cm[0][1]) *
        (cm[1][1] + cm[1][0]) *
        (cm[1][1] + cm[0][1])
    )
    return (cm[0][0]*cm[1][1] - cm[0][1]*cm[1][0])/d

# True positive rate
def get_tpr(cm):
    return float(cm[1][1]) / (cm[1][0] + cm[1][1])

# Positive predictive value (precision)
def get_ppv(cm):
    return float(cm[1][1]) / (cm[1][1] + cm[0][1])

def full_seq_computing(filename, th):
    cm = get_cm(filename, th, 2)
    print("USING E-VALUE OF THE FULL SEQUENCE...")
    q2 = get_q2(cm)
    print("tn=", cm[0][0], "fn=", cm[0][1])  # true negatives, false negatives
    print("fp=", cm[1][0], "tp=", cm[1][1])  # false positives, true positives
    mcc = get_mcc(cm)
    tpr = get_tpr(cm)
    ppv = get_ppv(cm)
    print("threshold=", th, "q2=", q2, "MCC=", mcc, "tpr=", tpr, "ppv=", ppv, "fullseq=", True)
    return

def single_domain_computing(filename, th):
    cm = get_cm(filename, th, 3)
    print("USING E-VALUE OF THE BEST DOMAIN...")
    q2 = get_q2(cm)
    print("tn=", cm[0][0], "fn=", cm[0][1])  # true negatives, false negatives
    print("fp=", cm[1][0], "tp=", cm[1][1])  # false positives, true positives
    mcc = get_mcc(cm)
    tpr = get_tpr(cm)
    ppv = get_ppv(cm)
    print("threshold=", th, "q2=", q2, "MCC=", mcc, "tpr=", tpr, "ppv=", ppv, "fullseq=", False)
    return

if __name__ == "__main__":
    filename = argv[1]
    th = float(argv[2])
    if len(argv) > 3:
        selection = int(argv[3])
        # selection = 1 → evaluate only full sequence
        # selection = 2 → evaluate only best domain
        # selection = 0 or undefined → evaluate both
    else:
        selection = 0
    if selection == 0:
        # Evaluate both full sequence and best domain using the same threshold
        full_seq_computing(filename, th)
        print("\n")
        single_domain_computing(filename, th)
    elif selection == 1:
        # Evaluate only full sequence
        full_seq_computing(filename, th)
    else:
        # Evaluate only best domain
        single_domain_computing(filename, th)