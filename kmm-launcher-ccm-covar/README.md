
A set of ec2 support scripts for the kmer-model which takes a bam and automatically generates fitted output to be placed into amazon's S3 storage system

## Prerequisites
Register an Amazon Elastic Cloud 2 (EC2) account following the instruction [here](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/get-set-up-for-amazon-ec2.html)


##	Quick example

Run the command in the git repo root:

```
docker pull thashim/kmm-launcher-cov
docker run --rm -w `pwd` -v /topfolder:/topfolder -i thashim/kmm-launcher-cov /kmm/run.onestrand.r example/covbinom.list example/auth.txt
```

There are three parameters changable in the running command:

+ `topfolder`: This should be the **common** top folder of the repo and all your bam/genome/list/auth files. For example if the repo and all you relevant files are under /cluster/project/wordfinder, you use either  */cluster* or */cluster/project* or */cluster/project/wordfinder*. The function of this argument is for the scripts inside Docker to access the data files on your file system.
+ `example/covbinom.list`: See step1
+ `example/auth.txt`: See step2

##Step1: Set up information needed to run on EC2

### auth.txt
Example: (/example/auth.txt)

```
realm:us-east-1
price:3.0
region:us-east-1d
rsa_key:/cluster/ec2/starcluster.rsa
access_key:REDACTED
secret_key:REDACTED
key_name:starcluster
mailaddr:thashim@csail.mit.edu
```


#####Useful options:

+ `price`: The max bid price. Setting this value to Inf will use on-demand allocation (cannot be killed) but will cost a fixed price of ~$.16 / hr. Use this setting only if there's heavy contention, and you cannot wait. $3 is reasonable. Set too low and your jobs will get killed before completed

+ `region`: The job submit regions. You can check the spot prices of a `c3.8xlarge` and pick a cheap region

+ `rsa_key`: The physical location of the key-pair file for your Amazon EC2 account. This file enables the user to remotely communicate with the EC2 instance created. Checkout [here](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html#having-ec2-create-your-key-pair) for a instruction to generate your key-pair file.

+ `access_key` and `secret_key`: The credentials for your Amazon EC2 account. Checkout [here](http://docs.aws.amazon.com/AWSSimpleQueueService/latest/SQSGettingStartedGuide/AWSCredentials.html) for instruction.

+ `mailaddr` sets the email address that gets emailed at the end of a job. The emails will probably get spam-boxed first, so check spam folder.

#####Optional options:

+ `itype`: The instance type: valid alternatives are cc2.8xlarge, to use this you must also change the AMI.

+ `ami`: The AMI type: you will want to use the HVM image (`ami-864d84ee`) if you use any other instances like cc2.8xlarge or r3.8xlarge.


##Step2: Set up parameters for the model
### *.list

Example: (/example/covbinom.list)

```
#bam.prefix /cluster/projects/wordfinder/bams/
#gbase /cluster/projects/wordfinder/data/genome/
#quality 0
#postfix .examplerun
#bucket_name batch_runs
#maxk 8
#k 1000
#resol 4
#read.max 5
#smooth.window 20
#genome mm10
#covariate dnase/bam1, dnase/bam2
#kbeta 100
#branch glm_v2
fos_run1,fos/bam1,fos/bam2,fos/bam3

#smooth.window 10
#quality 20
ctcf_run1,ctcf/bam1,ctcf/bam2,ctcf/bam3
```

Nearly all options can be overridden in a .list file.

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



#####Common arguments

+ `bam_prefix`: The top folder of **ALL** the bam files used, including ChIP-seq and covaraites DNase-seq bams.

+ `gbase`: The folder where genome files are stored. Do not change if run within gifford lab. Currently only hg19 and mm10 are supported, and their genome datafile can be found on [GERV website](http://gerv.csail.mit.edu).
+ `genome`: set to the organism genome. Currently only hg19 and mm10 are supported.

+ `bucket_name`: s3 bucket name. This should generally be your username / project name to avoid mixing multiple people's jobs. 

+ `postfix`: postfix applied to jobs. Each job will go into a S3 bucket where they are separated into folders named $experiment_name$$postfix$ (no spacer in between)


+ `branch`: Use *glm_v2* for the full model. Use *no91* for the model without DNase-seq covariates.

+ `covariate`: The path (relative to `bam_prefox`) of DNase-seq bams. 

+ `quality`: mapper quality cutoff, pick q=0 by default, q=20 if attempting to avoid repeat regions and other hard-to-map regions. q=0 was used in the GERV paper.

#####Tweakable parameters

+ `maxk`: Maximum kmer length to consider, 8 is generally good enough and the start of diminishing returns, and was used in the GERV paper.

+ `k`: the window size. The model looks within a `[-k,+k]` region around each Kmer match. Should be a multiple of RESOL. k=200 was used in the GERV paper.

+ `read.max`: truncate input at read.max to avoid giant read-spikes from affecting model. Generally 5-10 is good for experiments in the < 1 billion read range. 5 was used for the GERV paper.


+ `resol`: the resolution at which parameters are stored. for example, if K=1000, RESOL=5, then the model fits 200 paremters, each representing 5 bases. **RESOL MUST BE ABLE TO DIVIDE K**. resol=1 was used in the GERV paper

+ `smooth.window`: smooth the input by this many bases before feeding into the model. Useful for low-coverage experiments. Default of 10-20 is fine for all but extreme high or low coverage experiments. 50 was used in the GERV paper.

+ `kbeta`: window size for the prediction of target (ie ChIP) from covariate (ie DNase). **Don't include this parameter for model without DNase-seq covariates**. 200 was used in the GERV paper

#### Things to note:

+ **Make sure your specified bucket name exists in s3 before starting!**

+ **All the path of the bam files should be relative to the folder specificied after `bam_prefix`**

+ **Don't include `kbeta` and `covariates` for the model without DNase-seq covariates.**

