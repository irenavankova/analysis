export R=rg

export PREFIX=anvil:/lcrc/group/e3sm/ac.vankova/compass/sg_tests/sg_pull_w_fraz_yesC/
export SUFFIX=ocean/isomip_plus/planar/2km/z-star/Ocean0/simulation
export OUT=/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC

#declare -a hloc=("112" "122" "132" "142")
#declare -a sgr=("R" "A" "B" "C" "D" "E")
#declare -a hloc=("112")
#declare -a sgr=("N")
declare -a hloc=("112")
declare -a sgr=("R" "A" "B" "C" "N")

for h in "${hloc[@]}"
do
for s in "${sgr[@]}"
do

	export DIR=${R}/${R}_${h}${s}

	mkdir -p ${OUT}/${DIR}

	export IN=${PREFIX}/${DIR}/${SUFFIX}

	scp -r ${IN}/timeSeriesStatsMonthly.0002-12-01.nc ${OUT}/${DIR}
	scp -r ${IN}/restarts/restart.0003-01-01_00.00.00.nc ${OUT}/${DIR}
	scp -r ${IN}/init.nc ${OUT}/${DIR}

done
done

