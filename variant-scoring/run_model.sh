#!/bin/bash
docker pull haoyangz/wave-variant-scoring-covar;docker run --rm -w /scripts/ -v /cluster:/cluster -u $(id -u) -i haoyangz/wave-variant-scoring-covar /scripts/run.r /cluster/zeng/code/research/kmer/params_NFKB_covar preprocess > log/log_preprocess 2>&1

for i in `seq 1 22`;
do
	echo "docker pull haoyangz/wave-variant-scoring-covar;docker run --rm -w /scripts/ -v /cluster:/cluster -u $(id -u) -i haoyangz/wave-variant-scoring-covar /scripts/run.r /cluster/zeng/code/research/kmer/params_NFKB_covar score $i " | qsub -wd /cluster/zeng/code/research/kmer/ -o 'log/SGE_log/' -e 'log/SGE_log/' -N 'snp-scoring'
done

docker pull haoyangz/wave-variant-scoring;docker run --rm -w /scripts/ -v /cluster:/cluster -u $(id -u) -i haoyangz/wave-variant-scoring-covar /scripts/run.r /cluster/zeng/code/research/kmer/params_NFKB_covar combine.result > log/log_postprocess 2>&1 &
