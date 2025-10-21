mkdir -p /Users/peidli/work/data/SYMPHONY/SYM10RNAseq/fastq
cd /Users/peidli/work/data/SYMPHONY/SYM10RNAseq/
rsync -aP --append --partial \
  --files-from=fastq_list_subset.txt \
  sp41@farm22-head1.internal.sanger.ac.uk:/lustre/scratch125/casm/teams/team113/projects/imports/8094_Symphony101_Screen_Validation_Bulk_RNASeq/data/fastq/ \
  ./fastq/