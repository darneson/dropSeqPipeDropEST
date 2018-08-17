Description
------------------
This pipeline is based on [snakemake](https://snakemake.readthedocs.io/en/stable/) and the dropseq tools provided by the [McCarroll Lab](http://mccarrolllab.com/dropseq/). It starts from raw sequencing data and outputs digital gene expression matrices (DGEs).

This method is an adapted version of [dropSeqPipe](https://github.com/Hoohm/dropSeqPipe) which we have found to be much more sophisticated and flexible than our old [snakemake](https://snakemake.readthedocs.io/en/stable/) implementation of [Drop-seq Tools](http://mccarrolllab.com/download/1276/). The old version of the pipeline is available [here](https://github.com/darneson/DropSeq)

We have made a number of changes to the [dropSeqPipe](https://github.com/Hoohm/dropSeqPipe) workflow:

1). We have adapted the [dropSeqPipe](https://github.com/Hoohm/dropSeqPipe) workflow to work on a [SGE Cluster Environment](http://star.mit.edu/cluster/docs/0.93.3/guides/sge.html), specifically the [UCLA Hoffman2 Cluster](https://www.hoffman2.idre.ucla.edu/computing/sge/).

2). We have made the workflow compatible with the new pipeline [dropEST](https://github.com/hms-dbmi/dropEst) which is capable of handling reads aligning to intronic regions and results in a large boost in cell number and genes/UMIs per cell in our hands.

Updates
------------------
Uploaded 8-17-18
