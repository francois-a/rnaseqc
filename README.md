

## RNA-SeQC

**Note: this is an archival repository for v1.1.9 of the [RNA-SeQC](https://www.broadinstitute.org/cancer/cga/rna-seqc) quality control software for RNA sequencing data. For the latest version, please see https://github.com/broadinstitute/rnaseqc.**

#### Running the software

Java 1.7 is required.

RNA-SeQC provides two distinct functionalities:

1. Generation of quality control (QC) metrics.
2. Quantification of gene-level expression. This step should be performed using a collapsed transcript annotation (e.g., derived from GENCODE).

The commands for running RNA-SeQC in these two modes are:
```bash
# QC metrics
java -Xmx6g -jar RNA-SeQC.jar -n 1000 -s ${sample_id},${bam_file},${notes} -t ${annotation_gtf} -r ${genome_fasta} -o ${output_dir}

# Gene-level expression, based on a collapsed annotation GTF
java -Xmx6g -jar RNA-SeQC.jar -n 1000 -s ${sample_id},${bam_file},${notes} -t ${annotation_gtf} -r ${genome_fasta} -o ${output_dir} -noDoC -strictMode
```
For additional documentation, see http://www.broadinstitute.org/cancer/cga/rnaseqc_run
