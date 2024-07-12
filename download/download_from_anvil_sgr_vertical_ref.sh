
export PREFIX=blues:/lcrc/group/e3sm/ac.vankova/compass/sg_tests/sg_data_conserve_04_yesC/
export SUFFIX=ocean/isomip_plus/planar/2km/z-star/Ocean0/simulation
export OUT=/Users/irenavankova/Work/data_sim/SGR/idealized/sg_data_conserve_04_yesC

declare -a R=("f102" "f101" "f100")

for r in "${R[@]}"
do
	export DIR=${r}/${r}_10

	mkdir -p ${OUT}/${DIR}

	export IN=${PREFIX}/${DIR}/${SUFFIX}

	scp -r ${IN}/timeSeriesStatsMonthly.0002-12-01.nc ${OUT}/${DIR}
	scp -r ${IN}/restarts/restart.0003-01-01_00.00.00.nc ${OUT}/${DIR}
	scp -r ${IN}/init.nc ${OUT}/${DIR}

done

