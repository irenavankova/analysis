export R=f102

export PREFIX=blues:/lcrc/group/e3sm/ac.vankova/compass/sg_tests/sg_data_conserve_04_yesC/
export SUFFIX=ocean/isomip_plus/planar/2km/z-star/Ocean0/simulation
export OUT=/Users/irenavankova/Work/data_sim/SGR/idealized/sg_data_conserve_04_yesC

#declare -a hloc=("112" "122" "132" "142")
#declare -a sgr=("R" "A" "B" "C" "D" "E")
#declare -a hloc=("112")
#declare -a sgr=("N")
declare -a sgr=("" "A" "B" "C")
declare -a hloc=("110" "111" "112")

for h in "${hloc[@]}"
do
for s in "${sgr[@]}"
do

	export DIR=${R}/${R}_${h}000${s}

	mkdir -p ${OUT}/${DIR}

	export IN=${PREFIX}/${DIR}/${SUFFIX}

	scp -r ${IN}/timeSeriesStatsMonthly.0002-12-01.nc ${OUT}/${DIR}
	scp -r ${IN}/restarts/restart.0003-01-01_00.00.00.nc ${OUT}/${DIR}
	scp -r ${IN}/init.nc ${OUT}/${DIR}

done
done

