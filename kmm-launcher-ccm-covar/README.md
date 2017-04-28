Train a sequence-based model that predicts ChIP-seq or DNase-seq read counts from the sequence context.

## Prerequisites
+ [Docker](https://www.docker.com/)

##	Quick example

+ To run on Amazon EC2 cloud

	```
	docker pull haoyangz/gerv:launcher
	docker run --rm -w `pwd` -v TOPFOLDER:TOPFOLDER -i  \
			haoyangz/gerv:launcher /kmm/run.onestrand.r PARAM.LIST AUTH.TXT
	```

	+ `TOPFOLDER`: This should be the **common** top folder of the repo and all your bam/genome/list/auth files. For example if the repo and all you relevant files are under /cluster/project/wordfinder, you use either  */cluster* or */cluster/project* or */cluster/project/wordfinder*. The function of this argument is for the scripts inside Docker to access the data files on your file system.
    + `PARAM.LIST`: a file with running parameters. See step2 for details. 
    + `AUTH.TXT`: (optional for running on Amazon EC2 only) a file with EC2 account information. See step1 for details.

+ To run locally
	
	```
    docker pull haoyangz/gerv:launcher
    docker run --rm -w `pwd` -v TOPFOLDER:TOPFOLDER -v /etc/passwd:/etc/passwd -i -u $(id -u)  \
			haoyangz/gerv:launcher /kmm/run.onestrand.local.r PARAM.LIST
    ```

	+ `TOPFOLDER`: same as above
	+ `PARAM.LIST`: same as above

## Step1: Set up EC2-related configuration (Optional)
To run on Amazon EC2 cloud, set up the EC2-related configuration as instructed [here](https://github.com/gifford-lab/GERV/blob/master/kmm-launcher-ccm-covar/README_ec2.md).


## Step2: Set up parameters for the model

Examples: [local version](https://github.com/gifford-lab/GERV/blob/master/kmm-launcher-ccm-covar/example/param.local.list),[EC2 version](https://github.com/gifford-lab/GERV/blob/master/kmm-launcher-ccm-covar/example/param.list)


```
#bam.prefix /cluster/zeng/research/reads/
#gbase /cluster/projects/wordfinder/data/genome/
#quality 0
#bucket_name zengtest
#outdir /cluster/zeng/research/GERV_result/
#maxk 8
#k 100
#resol 4
#read.max 5
#genome hg19
#smooth.window 20

#postfix .withoutcov
CTCF_ENCODE_GM12878,CTCF/CTCF_crowford/ALL10/141104_Zeng_hg19_CTCF_GM12878_None_Rep1-ready.bam

#postfix .withcov
#kbeta 100
#covariate DNase/GM12878.DNase/141116_Crawford_hg19_DNase_GM12878_None_Rep1_bwa_hg19/141116_Crawford_hg19_DNase_GM12878_None_Rep1/141116_Crawford_hg19_DNase_GM12878_None_Rep1-ready.bam
CTCF_ENCODE_GM12878,CTCF/CTCF_crowford/ALL10/141104_Zeng_hg19_CTCF_GM12878_None_Rep1-ready.bam
```

The general format of a .list file is

```
#variable_name1 value
#variable_name2 value
[...]

experiment_name,bam_1.bam,bam_2.bam [..]

#variable_name1 value
[...]
experiment_name_2,bam_1 [...]
```

The launcher parses from top to bottom, setting each variable_name to value. Whenever it encounters a line without a `#` character, it will launch a KMM-job, assuming that the first entry is the experiment name and any following it are bams.

Later variable assignment lines starting with `#` will override earlier ones. In the example above, `fos_run1` launches with a `quality` parameter of 0, but `ctcf_1` is launched with `quality` of 20 due to the later override line.


##### Common arguments

+ `bam_prefix`: The top folder of **ALL** the bam files used, including ChIP-seq and covaraites DNase-seq bams.

+ `gbase`: The folder where genome files are stored. Do not change if run within gifford lab. Currently only hg19 and mm10 are supported, and their genome datafile can be found on [GERV website](http://gerv.csail.mit.edu).
+ `genome`: set to the organism genome. Currently only hg19 and mm10 are supported.

+ output related:
	+ Top directory
		+	(When running on EC2) `bucket_name`: s3 bucket name.  **Make sure your specified bucket name exists in s3 before starting!**. This should generally be your username / project name to avoid mixing multiple people's jobs. 

		+	(When running locally) `outdir`: The top directory to save all the input and outputfor each experiment. 

	+ Subdirectory
		+	`postfix`: a postfix appended to `experiment_name`. It's useful when you wish to try different sets of parameters on the same batch of jobs, where you can simply specify a  different `postfix` for each parameter set, for instance ".withcovariate" or ".defaultparams". Each job will go into a subfolder named `experiment_name`+`postfix` (+ denotes string concatenation).


+ `covariate`: The path (relative to `bam_prefox`) of DNase-seq bams. **Set to 'none' for model without DNase-seq covariates**

+ `quality`: mapper quality cutoff, pick q=0 by default, q=20 if attempting to avoid repeat regions and other hard-to-map regions. q=0 was used in the GERV paper.

##### Tweakable parameters

+ `maxk`: Maximum kmer length to consider, 8 is generally good enough and the start of diminishing returns, and was used in the GERV paper.

+ `k`: the window size. The model looks within a `[-k,+k]` region around each Kmer match. Should be a multiple of RESOL. k=200 was used in the GERV paper.

+ `read.max`: truncate input at read.max to avoid giant read-spikes from affecting model. Generally 5-10 is good for experiments in the < 1 billion read range. 5 was used for the GERV paper.


+ `resol`: the resolution at which parameters are stored. for example, if K=1000, RESOL=5, then the model fits 200 paremters, each representing 5 bases. **RESOL MUST BE ABLE TO DIVIDE K**. resol=1 was used in the GERV paper

+ `smooth.window`: smooth the input by this many bases before feeding into the model. Useful for low-coverage experiments. Default of 10-20 is fine for all but extreme high or low coverage experiments. 50 was used in the GERV paper.

+ `kbeta`: window size for the prediction of target (ie ChIP) from covariate (ie DNase). 200 was used in the GERV paper. This will be ignored if `covariate` is set as "none".

#### Things to note:

+ **All the path of the bam files, including the covariates, should be relative to `bam_prefix`**

## Understand the output
+	Log files
	+	`runlog.txt`: the run log for the training
	+	`nohup.txt`: the run log for every job, including preprocessing, training, postprocessing and so on.
	+	`runall.sh`: the command that the program runs.
	+	`param.list`: the parameter file you used
	
+	Input data

	All the input reads and accessory data are in 'input'

+	Trained model

	In `output` folder, the parameters for each training iteration (15 iter in total) are saved in binary format (4 byte double). In R you can read by `readBin(filename,double(),size=4,n=filesize)`
	+	`yall`: the influence profile of k-mers (size=(sum(4^(0:kmax)))*k2/resol, all the kmers concatenated). K-mers are ordered as "A,T,G,C,AA,AT,AG,AC,TA,TT,TG,TC..."
	+	`x0`: offsets (size=2)
	+	`eta`: regularization coefficient (size=1)
	+	`heldout`: loss evaluated on heldout data (size=1)

+	Evaluation

	In `summaries` folder, you can find predicted read counts and preliminaly quality analysis.
	+	`fitted.in`: the predicted read counts across genome using the parameters with the lowest heldout loss
	+	K-mer profile
		+	`profiles.pdf`: the influence profiles of k-mers
		+	`hctest-dist2.pdf`: profile comparision between the original sequence and the reverse-compliment of each k-mer
		+	`heatplot.pdf`: k-mer clustering based on influence profile.
	+	Perforamnce
		+	`heldout.pdf`: how the heldout loss changes with training iteration
		+	`cor.txt`: correlation
		+	`hits.txt`: coverage of hits
	+	Comparision between predicted vs. observed read count	
		+	`err.pdf`
		+	`detail.pdf`
		+	`callplots.pdf`

