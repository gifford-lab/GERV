The source code for GERV accompanying the following paper: 

Zeng, Haoyang, Hashimoto, Tatsunori, Kang, Daniel D & Gifford, David K (2015). GERV: A Statistical Method for Generative Evaluation of Regulatory Variants for Transcription Factor Binding. *Bioinformatics*. [[link]](http://bioinformatics.oxfordjournals.org/content/early/2015/10/15/bioinformatics.btv565.abstract)

##Pre-requisite
Packed with [Docker](https://www.docker.com/), GERV is universally runnable across different versions of linux systems without the need for further configuration. So the only pre-requisite is to install Docker.

## Running the model
The GERV model is consisted of two parts:

* Model training
* Variant Scoring


###Model training
In this part, GERV trains a k-mer model to fit the transcription factor binding ChIP-seq datasets that user provided. In the future, the trained model of a library of important transcription factors will be provided on the [GERV website](http://gerv.csail.mit.edu) for researchers who wish to skip this step and to score variants directly on existing models.

Please refer to the [README.md](https://github.com/gifford-lab/GERV/blob/master/kmm-launcher-ccm-covar/README.md) for more details on running this part.

###Variant Scoring
In this part, GERV use a k-mer model trained in the first step to score a set of variants that user provided.

Please refer to the [README.md](https://github.com/gifford-lab/GERV/blob/master/variant-scoring/README.md) for more details on running this part.




