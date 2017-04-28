# GERV
[![License](https://img.shields.io/hexpm/l/plug.svg)](LICENSE)

GERV is a statistical model that predicts the effect of genetic variants on transcription factor binding using ChIP-seq and DNase-seq data. 

## Pre-requisite
Packed with [Docker](https://www.docker.com/), GERV is universally runnable across different versions of linux systems without the need for further configuration. So the only pre-requisite is to install Docker.

## Running the model
The GERV model is consisted of two parts:

* Model training
* Variant Scoring


### Model training
In this part, GERV trains a k-mer model to fit the transcription factor binding ChIP-seq datasets that user provided. In the future, the trained model of a library of important transcription factors will be provided on the [GERV website](http://gerv.csail.mit.edu) for researchers who wish to skip this step and to score variants directly on existing models.

Please refer to the [README.md](https://github.com/gifford-lab/GERV/blob/master/kmm-launcher-ccm-covar/README.md) for more details on running this part.

### Variant Scoring
In this part, GERV use a k-mer model trained in the first step to score a set of variants that user provided.

Please refer to the [README.md](https://github.com/gifford-lab/GERV/blob/master/variant-scoring/README.md) for more details on running this part.


## License and Citation
GERV is released under [Apache License 2.0](https://github.com/gifford-lab/GERV/blob/master/LICENSE).

Please cite GERV if it helps in your research:

```
@article{Zeng2015,
author = {Zeng, Haoyang and Hashimoto, Tatsunori and Kang, Daniel D and Gifford, David K},
doi = {10.1093/bioinformatics/btv565},
journal = {Bioinformatics},
month = {oct},
title = {{GERV: A Statistical Method for Generative Evaluation of Regulatory Variants for Transcription Factor Binding}},
year = {2015}
}

```




