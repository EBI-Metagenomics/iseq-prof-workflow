# Baseline benchmark

- Download Pfam database.
- Download Pfam clans database.
- Download 980 bacteria and 20 archaea distinct whole-genome sequences
  from GenBank, randomly chosen.
- Run HMMER3 on coding sequences (CDSs) in the amino acid space
  (defined by GenBank). The high-quality annotations become the ground truth.
- Run iSeq on coding sequences in the DNA/RNA space (defined by GenBank).
  The resulting annotations are evaluated against the ground truth.
