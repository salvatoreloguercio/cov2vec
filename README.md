cov2vec is a systematic effort to obtain SARS CoV-2 genome embeddings by encoding viral genomes with protein language models - for specific applications (i.e. improve biomedically relevant ML tasks), and globally for learning meaningful representations of the viral genome as an evolving, high-dimensional genomic manifold.
Input: called mutations from a viral sequence deposited to GISAID (currently using gff3 files from CNCB; future support for outbreak.info / nextstrain input).
gff2fasta: converts called mutations back to a mutated sequence file in fasta format, by introducing amino-acid altering mutations (indels and missense mutations).
fasta2vec: uses a SOTA protein language model (ESM_1b or ESM_1v) pre-trained on a large protein sequence corpus (UniProt50, UniProt 90).
Output: per-protein and genome embedding for the input viral sequence.
