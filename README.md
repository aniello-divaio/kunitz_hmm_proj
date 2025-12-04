# In-Silico Modeling of the Kunitz-Type Domain Using Profile HMM
This repository presents a full workflow for building and evaluating a profile Hidden Markov Model (HMM) to identify the Kunitz protease inhibitor domain (PF00014) using structural and sequence data. The workflow integrates structural alignment, sequence filtering, model training, dataset construction, and performance evaluation.

---
## Tools and Scripts Used
- `cd-hit`, `hmmbuild` (HMMER), `BLAST+`, `hmmsearch`, `PDBeFold`
- `clean_fasta.sh`, `get_seq.py`, `performance.py`
- `performance_tsv.py`, `performance_graph.R`

## 1. Structural Data Acquisition and Preprocessing

### 1.1 Retrieve PDB entries with:
Go to (PDB)[https://www2.rcsb.org/] and with advanced research get a set of proteins with with Kunitz domain (N=160)
- PFAM ID = PF00014
- Resolution ≤ 3.5 Å
- Sequence length between 45 and 80 residues

The resulting sequences will be downloaded as a custom report (cvs file format) following these specifics:
- Entry ID, PDB ID, Resolution
- Sequence, Asym ID, PFAM annotation, Entity ID

`rcsb_pdb_custom_report.csv` will be converted in a FASTA file using the following bash command:
```
tr -d '"' < rcsb_pdb_custom_report.csv | awk -F ',' '{if (length($2)>0) {name=$2}; print name,$6,$8,$9}' | grep PF00014 | awk '{print ">"$1"_"$3; print $2}' > pdb_kunitz.fasta
```

### 1.2 Run `cd-hit`
Open the terminal and run 'cd-hit' with standard threshold  (sequence identity = 90%) allowed the identification of the set of proteins (N = 25):

```bash
cd-hit -i pdb_kunitz.fasta -o pdb_kunitz_cluster.txt    
```

**Results:**
- 160 sequences → 25 clusters

Two files were generated:
• A .clstr file (e.g. pdb_kunitz_cluster.txt.clstr)
• A .txt file (e.g. pdb_kunitz_cluster.txt — your clustered sequences)
The .clstr file contains raw cluster info, but its format may result hard to parse directly.

---

## 2. Multiple Structural Alignment
In order to perform MSA all the proteins IDs from the clustering were extracted, using the command:
```bash
grep '^>' pdb_kunitz_cluster.txt | sed 's/^>//' | sed 's/_/:/' > pdb_kunitz_ids_25.txt
```

### 2.1 Format input for PDBeFold:
[PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/) is the platform where the MSA was performed, with the following criteria:
- Flag Submission Form: multiple
- Upload the list `pdb_kunitz_ids_25.txt`

Submit your query and download the resulting `kunitz_msa_25.fasta`.
From this file three sequences were manually removed:
- 2ODY:F (highest number of aligned residue (Nres=127))
- 5JBT:Y (lowest number of aligned residues (Nres = 38))
- 4BQD:A removed due to the lack in the its chain A of one or both β-strands typical of Kunitz fold

Following the same procedure, perform a second multiple structural alignment using the set of 22 proteins, retrieving `kunitz_msa_22.fasta`. 

The downloaded FASTA file may contain inconsistencies such as, lowercase amino acids or extra information in the headers. To clean and standardize the file, use the `clean_fasta.sh` script provided in the repository:
``` bash
clean_fasta.sh kunitz_msa_22.fasta
```
`kunitz_msa_22.ali` will be used to build the HMM

## 3. Build the HMM Model
To build and run an HMM Model, go to your terminal and run: 

```bash
hmmbuild kunitz.hmm kunitz_MSA_clean.ali
```
Output file: `kunitz.hmm`

---

## 4. Validation Dataset Compilation
#### Collection of protein IDs to create the positive and the negative datasets:

### 4.1 Download the datasets:
From [UniProtKB](https://www.uniprot.org/) collect all human proteins containg a Kunitz domain from UniProtKB database (N = 18) with these filters: 
- Human (Taxonomy [OC] = 9606) 
- PFAM id = PF00014
- Reviewed: Yes

Output file: `kunitz_human.fasta`

Collect all not-human proteins containing a kunitz domain from the UniProtKB (N = 380) and download the fasta file (e.g kunitz_not_human.fasta) with these filters:
- Not human (NOT Taxonomy [OC] = 9606)
- PFAM id = PF00014
- SwissProt reviewed

 Download the fasta file `kunitz_not_human.fasta`. 

### 4.2 Merge datasets:

Merge the two retrived Kunitz domain dataset files to form a unified collection of positive examples to test the HMM, using the command:

```bash
cat kunitz_not_human.fasta kunitz_human.fasta > kunitz_all.fasta
```

---

## 5. Remove Highly Similar Sequences
Remove all sequences with a `sequence identity ≥ 95%` and a `Nres ≥ 50` mapping the aligned sequences set on the positive dataset. 

### 5.1 Build BLAST DB:
Create a blast database with the kunitz proteins from UniProtKB/SwissProt
```bash
makeblastdb -in kunitz_all.fasta -input_type fasta -dbtype prot -out kunitz_all.fasta
```

### 5.2 Run a blastp search on the aligned sequences:
```bash
blastp -query kunitz_MSA_clean.ali -db kunitz_all.fasta -out kunitz_pdb22.blast -outfmt 7
```
output file: `kunitz_pdb22.blast`

### 5.3 Filter highly similar hits and create a file with the ids to remove:
```
grep -v "^#" kunitz_pdb22.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u > high_match_22.txt
cut -d'|' -f2 high_match_22.txt > to_remove.txt
```

### 5.4 Extract the unmatched ids for the proteins that will form the positive database:
Create the file txt with all the ids
```bash
grep "^>" kunitz_all.fasta | cut -d'|' -f2 > kunitz_all.txt
comm -23 <(sort kunitz_all.txt) <(sort to_remove.txt) > kunitz_final.txt
```


## 6. Negative Dataset Preparation
Collect all SwissProt reviewed not-kunitz proteins (573.263) and download the fasta file (`uniprot_sprot.fasta`, the whole uniprot databses without the Kunitz proteins). Filters are:
- Not PFAM id PF00014
- SwissProt

### 6.1 Create a list of the IDs of UniProtKB/SwissProt dataset

```bash
grep ">" uniprot_sprot.fasta | cut -d "|" -f2 > uniprot_sprot.txt
```

### 6.2 Remove the positive IDs from the list of the negatives
```
comm -23 <(sort uniprot_sprot.txt) <(sort kunitz_final.txt) > negatives.txt
```
output file: `negatives.txt`

---

## 7. Dataset Subsetting
### 7.1 Random sorting of the ID files
```bash
sort -R kunitz_final.txt > kunitz_final_random.txt
sort -R negatives.txt > negatives_random.txt
```
### 7.2 Subsetting the randomized sets in halves
```

# Positives:
head -n 183 kunitz_final_random.txt > pos_1.txt
tail -n 182 kunitz_final_random.txt > pos_2.txt

# Negatives:
head -n 286416 negatives_random.txt > neg_1.txt
tail -n 286416 negatives_random.txt > neg_2.txt
```

### 7.3 Generate FASTA using `get_seq.py`:
```bash
# Positives:
python3 get_seq.py pos_1.txt kunitz_all.fasta > pos_1.fasta
python3 get_seq.py pos_2.txt kunitz_all.fasta > pos_2.fasta

# Negatives: 
python3 get_seq.py neg_1.txt uniprot_sprot.fasta > neg_1.fasta
python3 get_seq.py neg_2.txt uniprot_sprot.fasta > neg_2.fasta
```

---

## 8. HMM Scanning and Classification

Testing on positive dataset:
```bash
hmmsearch -Z 1000 --max --tblout pos_1.out kunitz.hmm pos_1.fasta
hmmsearch -Z 1000 --max --tblout pos_2.out kunitz.hmm pos_2.fasta

grep -v "^#" pos_1.out |awk '{split($1,a,"\|"); print a[2],1,$5,$8}' |tr " " "\t" >pos_1.class
grep -v "^#" pos_2.out |awk '{split($1,a,"\|"); print a[2],1,$5,$8}' |tr " " "\t" >pos_2.class
```
Testing on negative dataset:
```
hmmsearch -Z 1000 --max --tblout neg_1.out kunitz.hmm neg_1.fasta
hmmsearch -Z 1000 --max --tblout neg_2.out kunitz.hmm neg_2.fasta

grep -v "^#" neg_1.out |awk '{split($1,a,"\|"); print a[2],0,$5,$8}' |tr " " "\t" >neg_1.class
grep -v "^#" neg_2.out |awk '{split($1,a,"\|"); print a[2],0,$5,$8}' |tr " " "\t" >neg_2.class

comm -23 <(sort neg_1.txt) <(cut -f 1 neg_1.class | sort) | awk '{print $1"\t0\t10.0\t10.0"}' >> neg_1.class
comm -23 <(sort neg_2.txt) <(cut -f 1 neg_2.class | sort) | awk '{print $1"\t0\t10.0\t10.0"}' >> neg_2.class
```
Combine the two datasets:
```
cat pos_1.class neg_1.class >set_1.class
cat pos_2.class neg_2.class >set_2.class
```
---

## 9. Performance Evaluation

Use `performance.py`:
```bash
for i in `seq 1 12`; do python3 performance.py set_2.class 1e-$i; done
for i in `seq 1 12`; do python3 performance.py set_1.class 1e-$i; done
```
### 9.1 Visualize the Performance
In order to visualize how HMM performed, `performance_tsv.py` was used, following these commands:
```
python3 performance_tsv.py set_1.class
python3 performance_tsv.py set_2.class
```
After obtained .tsv files, those were plotted using Rstudio, specifically the script `perfomance_graph.R`
### Resulting Plot
![Plot](9_Performance_Evaluation/Rplot.png)

### Confusion Matrix
`set_1.class`
Best-Domain (E-value = 1e-06)
| TN         | TP      | FP    | FN    |
| ---------- | ------- | ----- | ----- |
| **286416** | **182** | **0** | **1** |

| Metric                | Value     |
| --------------------- | --------- |
| **MCC**               | 0.9972623 |
| **Q2 (Accuracy)**     | 0.9999965 |
| **TPR (Sensitivity)** | 0.9945355 |
| **TNR (Specificity)** | 1.0       |
| **PPV (Precision)**   | 1.0       |
| **FNR**               | 0.0054645 |
| **FPR**               | 0.0       |



`set_2.class`
Best-Domain (E-value = 1e-05)
| TN         | TP      | FP    | FN    |
| ---------- | ------- | ----- | ----- |
| **286416** | **179** | **0** | **3** |

| Metric                | Value     |
| --------------------- | --------- |
| **MCC**               | 0.9917188 |
| **Q2 (Accuracy)**     | 0.9999895 |
| **TPR (Sensitivity)** | 0.9835165 |
| **TNR (Specificity)** | 1.0       |
| **PPV (Precision)**   | 1.0       |
| **FNR**               | 0.0164835 |
| **FPR**               | 0.0       |
