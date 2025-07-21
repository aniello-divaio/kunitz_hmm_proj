# Analisi HMM per il dominio Kunitz

Questo progetto costruisce un modello HMM a partire da un allineamento strutturale.

## Requisiti

- `cd-hit`
- `hmmbuild`
- `hmmsearch`

## Comandi

```bash
cd-hit -i input.fasta -o clustered.fasta -c 0.9
hmmbuild model.hmm alignment.fasta
hmmsearch model.hmm test_set.fasta

