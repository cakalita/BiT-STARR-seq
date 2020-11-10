# BiT-STARR-seq
code here is related to work from "High-throughput characterization of genetic effects on DNA–protein binding and gene transcription" https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6211638/

alignandprep.sh takes fastq through alignment and removing duplicates via UMI tools, to getting allelic counts. The final counts file dedup.${EXP}.q20.query.txt then goes into R for analysis.

bitstarr_base_script.R contains script to process DNA library and then individual ASE and meta-analysis ASE.
