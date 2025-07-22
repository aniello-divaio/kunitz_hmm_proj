# In Silico Modeling of the Kunitz-Type Domain Using Profile Hidden Markov Models
#### This repository presents a full workflow for building and evaluating a profile Hidden Markov Model (HMM) to identify the Kunitz protease inhibitor domain (PF00014) using structural and sequence data. The workflow integrates structural alignment, sequence filtering, model training, dataset construction, and performance evaluation.

---

## 1. Structural Data Acquisition and Preprocessing

### 1.1 Retrieve PDB entries with:
Go to https://www2.rcsb.org/ and with advanced research get a set of proteins with with Kunitz domain (N=160)
- Resolution ≤ 3.5 Å
- PFAM ID = PF00014
- Sequence length between 45 and 80 residues

### 1.2 Create a custom PDB report:
Following these specifics:
- Entry ID, PDB ID, Resolution
- Sequence, Asym ID, PFAM annotation, Entity ID
Download it in a  `.cvs` file, then convert it in FASTA using `convert_to_fasta.sh`.

### 1.3 Run `cd-hit`
Open the terminal and run 'cd-hit' with 90% threshold:

```bash
cd-hit -i kunitz_sequences.fasta -o kunitz_clustered.fasta -c 0.9
```

**Results:**
- 160 sequences → 25 clusters

---

### 1.4 Extract representative sequences:
Extract the representative sequence IDs form the '.clstr' file with:

```bash
grep '*' kunitz_clustered.fasta.clstr | cut -d '>' -f2 | cut -d '.' -f1 > representative_ids.txt
```
Retrieve the corresponding sequences from the original FASTA file and save them to a new file:
```
for i in $(cat representative_ids.txt); do
  grep -A 1 "^>$i" kunitz_sequences.fasta | tail -n 2 >> kunitz_representatives.fasta
done
```

---

## 2. Multiple Structural Alignment

### 2.1 Format input for PDBeFold:
Modify the previous file with: 
```bash
grep '^>' pdb_kunitz_cluster.txt | sed 's/^>//' | sed 's/_/:/' > pdb_kunitz_ids_25.txt
```
to obtain a file with only codes for the protein to input inside the alignment tool. 

Go to [PDBeFold](https://www.ebi.ac.uk/msd-srv/ssm/) -> Launch PDBeFold -> Flag Submission Form: multiple -> Upload the list `pdb_kunitz_ids_25.txt` -> Submit your query. 
Download `efold_output.txt`, then clean it using the script:

```bash
./clean_fasta.sh
```

Output: `kunitz_hmm_ready.fasta`

---

## 3. Build the HMM Model
To build and run an HMM Model, go to your terminal and run: 

```bash
hmmbuild kunitz_domain.hmm kunitz_hmm_ready.fasta
```

---

## 4. Validation Dataset Compilation

### 4.1 Download:
- Human Kunitz proteins (N = 18)
- Non-human Kunitz proteins (N = 380)
- Both reviewed, PF00014

### 4.2 Merge datasets:
```bash
cat NOT_human_kunitz.fasta kunitz_all.fasta > hmm_test_set.fasta
```

---

## 5. 🚫 Remove Highly Similar Sequences

### 5.1 Build BLAST DB:
```bash
makeblastdb -in hmm_set.fasta -dbtype prot
```

### 5.2 Run BLAST and filter:
```bash
blastp -query kunitz_hmm_ready.ali -db hmm_set.fasta -outfmt 7 > kunitz_pdb22.blast

grep -v "^#" kunitz_pdb22.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u > high_match_22.txt
cut -d'|' -f2 high_match_22.txt > to_remove.txt
```

### 5.3 Extract final positive IDs:
```bash
grep "^>" hmm_set.fasta | cut -d'|' -f2 > kunitz_all.txt
comm -23 <(sort kunitz_all.txt) <(sort to_remove.txt) > kunitz_final.txt
```

---

## 6. 📁 Negative Dataset Preparation

Download SwissProt without PF00014:
```bash
grep ">" uniprot_sprot.fasta | cut -d "|" -f2 > uniprot_sprot.txt
comm -23 <(sort uniprot_sprot.txt) <(sort kunitz_final.txt) > negatives.txt
```

---

## 7. 🧪 Dataset Subsetting

```bash
sort -R kunitz_final.txt > kunitz_final_random.txt
sort -R negatives.txt > negatives_random.txt

# Split sets
head -n 183 kunitz_final_random.txt > pos_1.txt
tail -n 182 kunitz_final_random.txt > pos_2.txt

head -n 286416 negatives_random.txt > neg_1.txt
tail -n 286416 negatives_random.txt > neg_2.txt
```

### Generate FASTA using `get_seq.py`:
```bash
python3 get_seq.py pos_1.txt NOT_uniprot_sprot.fasta > pos_1.fasta
python3 get_seq.py pos_2.txt NOT_uniprot_sprot.fasta > pos_2.fasta

python3 get_seq.py neg_1.txt uniprot_sprot.fasta > neg_1.fasta
python3 get_seq.py neg_2.txt uniprot_sprot.fasta > neg_2.fasta
```

---

## 8. 🔍 HMM Scanning and Classification

```bash
hmmsearch -Z 1000 --max --tblout pos_1.out kunitz_domain.hmm pos_1.fasta
hmmsearch -Z 1000 --max --tblout neg_1.out kunitz_domain.hmm neg_1.fasta

# Format
grep -v "^#" pos_1.out | awk '{split($1,a,"|"); print a[2],1,$5,$8}' | tr " " "\t" > pos_1.class
grep -v "^#" neg_1.out | awk '{split($1,a,"|"); print a[2],0,$5,$8}' | tr " " "\t" > neg_1.class

# Add missing negatives
comm -23 <(sort neg_1.txt) <(cut -f 1 neg_1.class | sort) | awk '{print $1"\t0\t10.0\t10.0"}' >> neg_1.class

# Combine sets
cat pos_1.class neg_1.class > set_1.class
```

---

## 9. 📊 Performance Evaluation

Use `performance.py`:
```bash
python3 performance.py set_1.class 1e-1
```

Loop version:
```bash
for i in $(seq 1 12); do
  echo "Threshold: 1e-$i"
  python3 performance.py set_1.class 1e-$i
done
```

Metrics reported:
- TP, FP, TN, FN
- MCC, Q2, TPR, PPV

---

## 🧪 Tools and Scripts Used

- `cd-hit`, `hmmbuild` (HMMER), `BLAST+`, `hmmsearch`, `PDBeFold`
- `convert_to_fasta.sh`, `clean_fasta.sh`, `get_seq.py`, `performance.py`
