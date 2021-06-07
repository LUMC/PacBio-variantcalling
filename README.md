# PacBio-variantcalling

![Tests](https://github.com/LUMC/PacBio-variantcalling/workflows/Continuous%20Integration/badge.svg)
[![image](https://img.shields.io/github/release/LUMC/PacBio-variantcalling.svg)](https://github.com/LUMC/PacBio-variantcalling/releases)
[![image](https://img.shields.io/github/release-date/LUMC/PacBio-variantcalling.svg)](https://github.com/LUMC/PacBio-variantcalling/releases)
[![DOI](https://zenodo.org/badge/314152152.svg)](https://zenodo.org/badge/latestdoi/314152152)

------------------------------------------------------------------------

## Documentation
To download the pipeline and all associated files, you can run
```bash
git clone https://github.com/LUMC/PacBio-variantcalling.git
cd PacBio-variantcalling
git submodule update --init --recursive
```

Next, install the
[conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
that is used to execute the pipeline, and activate the environment
```bash
conda env create --file environment.yml
conda activate PacBio-variantcalling
```

You can test the workflow using the following two commands. The first command
runs the sanity checks to make sure everything required is installed. The
second command will run the integration tests.

```bash
pytest --kwd tests --tag sanity
pytest --kwd tests --tag integration
```

To get a better sense of the pipeline and its inputs, you can manually run the
most general test case
```bash
cromwell run \
    --options tests/data/config/cromwell.options.json \
    --inputs tests/data/config/variant_calling.json \
    PacBio-variantcalling.wdl
```

This will run the pipeline for you using the
`tests/data/config/variant_calling.json` example configuration file. After the
pipeline has completed, you can find the full execution folder in
`cromwell-executions`. The workflow outputs have also
been copied to the `test-output` folder in the current directory, as is
specified in the `tests/data/config/cromwell.options.json` options file.

### Pipeline configuration file
To generate an input configuration file for the PacBio pipeline, please run the
following command.

```bash
womtool inputs --optional-inputs false PacBio-variantcalling.wdl
{
  "VariantCalling.barcodesFasta": "File",
  "VariantCalling.referenceFileDict": "File",
  "VariantCalling.referenceFileIndex": "File",
  "VariantCalling.referenceFile": "File",
  "VariantCalling.subreadsFile": "File",
  "VariantCalling.referencePrefix": "String"
}
```

If you also want to see the optional pipeline inputs, you can leave out the
`--optional-inputs false` argument.

### Common configuration options
| Setting                           | Type | Required | Description |
| --------------------------------- | ---- | -------- | ----------- |
| VariantCalling.barcodesFasta      | File | Required | Fasta file with the barcodes to be used by Lima for demultiplexing. |
| VariantCalling.referenceFileDict  | File | Required | The picard dictionary file for the reference. |
| VariantCalling.referenceFileIndex | File | Required | The samtools index file for the reference. |
| VariantCalling.referenceFile      | File | Required | The fasta reference file. |
| VariantCalling.subreadsFile       | File | Required | The subreads bam file. |
| VariantCalling.referencePrefix    | String | Required | The name of the reference. |
| VariantCalling.useDeepVariant     | Boolean | Optional | Use DeepVariant instead of GATK4 for variant calling. |
| VariantCalling.generateGVCF       | Boolean | Optional | Generate g.vcf files for all sample. This is extremely slow when used in combination with `VariantCalling.useDeepVariant`. |
| VariantCalling.subreadsIndexFile  | File | Optional | The index for the subreads bam file. If not specified, the pipeline will index the subreads file, this takes approximately two hours on 600GB of data. |
| VariantCalling.limaBarcodes       | File | Optional | TSV file containing the mapping from barcodes to sample names (forward_barcode reverse_barcode sample_name). This is used by MultiQC to rename the barcodes to their apropriate sample names. |
| VariantCalling.targetGenes        | File | Optional | Bed file containing the target genes. Used to determine the PGx phasing and Picard HsMetrics. |
| VariantCalling.dbsnp              | File | Optional | dbSNP file used to annotate the discovered variants. The results are displayed in the MultiQC report. |
| VariantCalling.dbsnpIndex         | File | Optional | Index for the dbSNP file, required when VariantCalling.dbsnp is specified. |

### Check your configuration file
If you have create your own configuration file, you can use the following
command to make sure all inputs are valid. Replace
`tests/data/config/variant_calling.json` with the path to your own
configuration file.

```bash
womtool validate --inputs tests/data/config/variant_calling.json PacBio-variantcalling.wdl
Success!
```
