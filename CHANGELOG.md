Changelog
==========

<!--
Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

v2.1.0-dev
---------------------------
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
