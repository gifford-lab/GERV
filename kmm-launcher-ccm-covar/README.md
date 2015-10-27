
A set of ec2 support scripts for the kmer-model which takes a bam and automatically generates fitted output to be placed into amazon's S3 storage system

**Make sure your specified bucket name exists in s3 before starting!**


##Example

Run the command in the git repo root:

```
docker pull thashim/kmm-launcher-cov/
docker run --rm -w `pwd` -v /cluster:/cluster -i thashim/kmm-launcher-cov/ /kmm/run.onestrand.r example/covbinom.list /cluster/ec2/auth.txt
```

##Configuring the KMM

### auth.txt
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
This file specifies information related to launch computing instance on Amazon EC2 cloud.

#####Useful options:

`price` sets the max bid price. Setting this value to Inf will use on-demand allocation (cannot be killed) but will cost a fixed price of ~$.16 / hr. Use this setting only if there's heavy contention, and you cannot wait. $3 is reasonable. Set too low and your jobs will get killed before completed

`region` sets the job submit regions, you can check the spot prices of a `c3.8xlarge` and pick a cheap region

`rsa_key` points to the physical location of the key-pair file for your Amazon EC2 account. This file enables the user to remotely communicate with the EC2 instance created. Checkout [here](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html#having-ec2-create-your-key-pair) for a instruction to generate your key-pair file.

`access_key` and `secret_key` are the credentials for your Amazon EC2 account. Checkout [here](http://docs.aws.amazon.com/AWSSimpleQueueService/latest/SQSGettingStartedGuide/AWSCredentials.html) for instruction.

`mailaddr` sets the email address that gets emailed at the end of a job. The emails will probably get spam-boxed first, so check spam folder.

#####Optional options:

`itype` sets the instance type: valid alternatives are cc2.8xlarge, to use this you must also change the AMI.

`ami` sets the AMI type: you will want to use the HVM image (`ami-864d84ee`) if you use any other instances like cc2.8xlarge or r3.8xlarge.

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

`bam_prefix`: the path to where bam files are stored. Launcher will look for `bam_prefix+bam_name` where `bam_name` is the name in the `experiment_name`.

`gbase`: path to where genome files are stored. do not change if run within gifford lab.

`quality`: mapper quality cutoff, pick q=0 by default, q=20 if attempting to avoid repeat regions and other hard-to-map regions.

`postfix`: postfix applied to jobs. Each job will go into a S3 bucket where they are separated into folders named `experiment_name+postfix`

`bucket_name`: s3 bucket name. This should generally be your username / project name to avoid mixing multiple people's jobs.

`genome`: set to the organism genome. Currently only hg19 and mm10 are supported.

`branch`: Use `glm_v2` for the full model. Use `no91` for the model without DNase-seq covariates.

`covariate`: list of experiments to use to predict the target experiment (ie DNase-seq predicting a ChIP). **Don't include this parameter for the model without DNase-seq covariates.**

#####Tweakable parameters

`maxk`: Maximum kmer length to consider, 8 is generally good enough and the start of diminishing returns.

`k`: the window size. The model looks within a `[-k,+k]` region around each Kmer match. Should be a multiple of RESOL.

`read.max`: truncate input at read.max to avoid giant read-spikes from affecting model. Generally 5-10 is good for experiments in the < 1 billion read range


`resol`: the resolution at which parameters are stored. for example, if K=1000, RESOL=5, then the model fits 200 paremters, each representing 5 bases. **RESOL MUST BE ABLE TO DIVIDE K**

`smooth.window`: smooth the input by this many bases before feeding into the model. Useful for low-coverage experiments. Default of 10-20 is fine for all but extreme high or low coverage experiments.

`kbeta`: window sie for the prediction of target (ie ChIP) from covariate (ie DNase). **Don't include this parameter for model without DNase-seq covariates.**

