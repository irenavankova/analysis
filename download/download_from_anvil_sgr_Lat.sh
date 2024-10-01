export R=rd

export FLDRSFX=L85
export FLDR=${R}${FLDRSFX}


export PREFIX=anvil:/lcrc/group/e3sm/ac.vankova/compass/sg_tests/sg_pull_w_fraz_yesC/${FLDR}
export SUFFIX=ocean/isomip_plus/planar/2km/z-star/Ocean0/simulation
export OUT=/Users/irenavankova/Work/data_sim/SGR/idealized/sg_pull_w_fraz_yesC/${FLDR}/${R}

#declare -a hloc=("112" "122" "132" "142")
#declare -a sgr=("R" "A" "B" "C" "D" "E")
#declare -a hloc=("112")
#declare -a sgr=("N")
declare -a hloc=("122")
declare -a sgr=("D" "A" "B" "C" "N")
#declare -a sgr=("N")

for h in "${hloc[@]}"
do
for s in "${sgr[@]}"
do

	export DIR=${R}_${h}${s}_${FLDRSFX}

	mkdir -p ${OUT}/${DIR}

	export IN=${PREFIX}/${DIR}/${SUFFIX}

	scp -r ${IN}/timeSeriesStatsMonthly.0002-12-01.nc ${OUT}/${DIR}
	scp -r ${IN}/restarts/restart.0003-01-01_00.00.00.nc ${OUT}/${DIR}
	scp -r ${IN}/init.nc ${OUT}/${DIR}

done
done

