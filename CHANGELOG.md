Changelog
==========

<!--
Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

v2.0.0-dev
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
