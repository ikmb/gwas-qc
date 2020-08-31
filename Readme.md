# GWAS Quality Control Pipeline User Manual

## Prerequisites
- Nextflow: https://www.nextflow.io/
- Java 8 or higher
- Singularity 3.4 or higher
- A to-be-QCed dataset in Plink format (bim/bed/fam). You can use the included example for test runs.

Please ensure that you have 16 GB RAM installed on the computer where you intend to run the pipeline (i.e. your local computer or your HPC compute nodes).

## Quick Start

1. Get the example dataset: https://raw.githubusercontent.com/ikmb/gwas-qc/master/Example.tar.gz
2. Unpack it in your home folder: ```tar xvaf Example.tar.gz -C $HOME; cd $HOME/Example```
3. Launch the pipeline: ```nextflow run -c pipeline.config ikmb/gwas-qc```

The pipeline output and reports will be written to the ```output``` directory.

## How to use

In most cases, the example configuration consisting of three files, ```QC.config```, ```dataset.config``` and ```pipeline.config``` can be used for your own dataset with only minor changes. To prepare your own dataset, perform the following steps:

1. Create a new directory to host your configuration and results
2. Copy the aforementioned configuration files to your newly-created directory
3. Adjust the names and locations of datasets and config files in your ```pipeline.config```. Note that absolute path names are required, relative names are not supported.
4. Change the items ```collection_name``` and ```allowed_diagnoses``` in your dataset.config according to your data. The collection name is the name prefix that will be used for output files and is arbitrary, as long as it would be a valid filename. The ```allowed_diagnoses``` key is used to specify which samples to filter based on specified diseases in case you have a multi-disease dataset. Disease names are arbitrary but should not contain whitespaces.
5. Create individual annotations (see below) for your dataset and place it among your source files. If your source files are named ```MyGWAS.bim/bed/fam```, the annotations are expected to be in ```MyGWAS_individuals_annotation.txt```.
6. Run the pipeline with ```nextflow run -c pipeline.config ikmb/gwas-qc```

### Individual Annotations

For proper analysis of batch and principal component effects, the pipeline makes use of certain sample information that you need to supply in a text file. It contains a single header line and one line for each sample in your dataset. The columns are separated by one or more whitespaces (or tab characters). The format is as follows:
```
familyID    individualID	paternalID	maternalID	sex	phenotype	batch	ethnicity_predicted	diagnosis	country
HG00096	    HG00096	        0	        0               2       2	        1000G	European	        Control	        Somewhere
HG00097	    HG00097	        0	        0	        1	2	        1000G	European	        Control	        Somewhere
HG00099	    HG00099	        0	        0	        2	2	        1000G	European	        Control	        Somewhere
HG00100	    HG00100	        0	        0       	1	1	        1000G	European	        Case	        Somewhere
HG00101	    HG00101	        0	        0	        2	1	        1000G	European	        Case	        Somewhere
HG00102	    HG00102	        0	        0	        2	2	        1000G	European	        Control	        Somewhere
```

* familyID, individualID, paternalID, maternalID: can be copied from the Plink FAM file. The ```individualID``` must be unique. ```paternalID``` and ```maternalID``` are required but are currently not used.
* sex, phenotype: same encoding as in the FAM file. For sex, 1/2 is male/female and for phenotype, 1/2 is control/case.
* batch: used by principal component analysis to find batch effects. Can be set to the same value for all samples if you only have one batch
* ethnicity_predicted: used to provide a reference frame for PCA plots. Currently, only ```European``` is supported
* diagnosis: if the ```phenotype``` is 1, this should be control. Otherwise, pick a disease name for the cases. This is used in PCA effect analysis and diagnosis filtering (see "How to use")
* country: the probably self-reported origin of the sample. Is used in PCA to show batch effects

Note that empty lines or comment lines are currently not supported.

## Advanced Configuration

### Local and Side-wide Configuration

You can extend the Nextflow default configuration by using (or creating) the following configuration files:
- a ```nextflow.config``` in the current project folder. This file is automatically read if Nextflow is launched from this folder.
- a site-wide ```$HOME/.nextflow/config```. This file is automatically read every time you run Nextflow

It is usually a good idea to place additional configuration items to the side-wide configuration (see below for examples).

### Shared Singularity Cache

In a shared compute enviroment such as HPC clusters, it is often useful to share the singularity cache with other Nextflow/Singularity users so they would not have to download the same large files over and over again. To specify the singularity cache directory, add the following line to your config:
```
singularity.cacheDir = '/some/shared/place/singularity'
```
Note that the directory must be accessible from all compute nodes.

### HPC Resources and Job Queues

By default, all processes are launched on the computer where the QC is started. This is usually not appropriate on HPC login nodes where jobs should be sheduled on different nodes. Nextflow provides support for a broad range of job submission systems, such as SLURM, SGE or PBS/Torque. Please review the [Nextflow documenation on compute resources](https://www.nextflow.io/docs/latest/executor.html).

For example, if you intend to use a SLURM environment, place the following code in your config (see "Shared Singularity Cache"):
```
process.executor = "slurm"
executor.queue = "partition-name"  // optional
executor.queueSize = 150           // optional, override max queue length for Nextflow
process.clusterOptions = '--qos=someval' // optional, if you need to supply additional args to sbatch
```

### Mounting Paths into the Singularity Container

To separate the operating system within the singularity container from the host system, Singularity only makes your home folder accessible to the insides. If your data files are stored in a directory different from your home (e.g. a shared storage ```/data_storage```), you will need to explicitly make it accessible by specifying additional singularity options in your config file:
```
singularity.runOptions = "-B /data_storage -B /some_other_storage -B /even_more_storage"
```

