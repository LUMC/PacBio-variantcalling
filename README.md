# PacBio-variantcalling
Variant calling workflow for targeted PacBio read data.

## Documentation
To download the pipeline and all associated files, you can run
```bash
git clone git@github.com:biowdl/PacBio-variantcalling.git
git submodule update --init --recursive
```

You can install the
[conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
that is used to execute the pipeline using
```bash
cd PacBio-variantcalling
conda install --file requirements-test.txt
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

For your own analysis, you can take the `variant_calling.json` file as a
starting point, and update the `subreadsFile` and `barcodesFasta` to point to
your own data. It is probably also best to leave out
`"VariantCalling.ccsChunks": "1"`, since that prevents running the CCS step in
parallel. By default, the pipeline will split the CSS step in 20 chunks to
speed up the analysis.
