Changelog
==========

<!--
Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

v3.0.0
---------------------------
+ CSS calling and demultiplexing with Lima have been removed from the pipeline.
  PacBio-variantcalling now expects one or more post-CCS bam files per sample,
  and only performs variantcalling.

v2.1.2
---------------------------
+ Allow nested inputs for better control over individual jobs.

v2.1.1
---------------------------
+ Fix bug in MultiQC pbmarkdup module for long library names.
+ Add VCF Metrics to MultiQC report (only used when VariantCalling.dbsnp is
specified).

v2.1.0
---------------------------
+ Include pbmarkdup in pipeline.
+ Upgrade MultiQC_PGx to latest version.
+ Update documentation and add conda environment file.

v2.0.2
---------------------------
+ Dummy release for Zenodo.

v2.0.1
---------------------------
+ Move pipeline out of BioWDL.
+ Make g.vcf output optional, because generating a g.vcf with DeepVariant is
extremely slow. This changes the default behaviour.

v1.0.0
---------------------------
+ Set scatter size for GATK and DeepVariant to 350 million basepairs.
+ Add scattering to the genotyping steps for GATK and DeepVariant.
+ Automatically convert target and bait bed files to picard .interval files,
    simplifying the output.
+ Add custom MultiQC configuration file.
+ Update tasks and the input/output names.
