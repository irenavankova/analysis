export R=FS

export PREFIX=anvil:/lcrc/group/e3sm/ac.vankova/compass/tests/fw_tracers/sgr_tracers_mpaso/
export SUFFIX=ocean/isomip_plus/planar/2km/z-star/Ocean0/simulation
export OUT=/Users/irenavankova/Work/data_sim/SGR/idealized/sgr_tracers_mpaso


declare -a sgr=("011" "001" "010" "110" "111")


for s in "${sgr[@]}"
do
	export DIR=${R}_${s}

	mkdir -p ${OUT}/${DIR}

	export IN=${PREFIX}/${DIR}/${SUFFIX}

	scp -r ${IN}/timeSeriesStatsMonthly.0001-12-01.nc ${OUT}/${DIR}
	scp -r ${IN}/restarts/restart.0002-01-01_00.00.00.nc ${OUT}/${DIR}
	scp -r ${IN}/init.nc ${OUT}/${DIR}

done

