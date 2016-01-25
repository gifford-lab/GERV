#!/bin/bash
docker pull haoyangz/gerv;docker run --rm -w /scripts/ -v /cluster:/cluster -u $(id -u) -i haoyangz/gerv /scripts/run.r /cluster/mypath/params_NFKB_covar preprocess > log/log_preprocess 2>&1

for i in `seq 1 22`;
do
	echo "docker pull haoyangz/gerv;docker run --rm -w /scripts/ -v /cluster:/cluster -u $(id -u) -i haoyangz/gerv /scripts/run.r /cluster/mypath/params_NFKB_covar score $i " | qsub
done

docker pull haoyangz/gerv;docker run --rm -w /scripts/ -v /cluster:/cluster -u $(id -u) -i haoyangz/gerv /scripts/run.r /cluster/mypath/params_NFKB_covar combine.result > log/log_postprocess 2>&1 &
