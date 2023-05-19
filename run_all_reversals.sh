#!/bin/bash

declare lat_min
declare lat_max
declare original_lat
declare time
declare dummy_input
declare executable
declare inputfile_root
declare infile
declare outfile
declare outfile_prefix
declare -i njob=0
declare -i nmax=27
declare -i njobs=1

lat_min=-85
lat_max=85
delta=5
inputfile_root="TAPE1"
original_lat="XYZAB"
dummy_input="TAPE1_dummy_lat_5x5_vertical"
executable="tji95_reversal.out"
magneticfield="reversal"
ddir="vert"
#declare -a time_arr=("766.958496" "771.506653" "773.321350" "774.515259" "775.504517" "779.465942" "782.647339" "783.786682" "787.197754")
declare -a time_arr=("794.474792" "794.247437" "794.020020" "793.792603" "793.565186" "793.337769" "793.110352" 
                     "792.882935" "792.655579" "792.428162" "792.200745" "791.973328" "791.745911" "791.518494" 
		     "791.291138" "791.063721" "790.836304" "790.608887" "790.381470" "790.379211" "790.154053" 
		     "789.926697" "789.699280" "789.471863" "789.244446" "789.017029" "788.789612" "788.562195" 
		     "788.334839" "788.107422" "787.880005" "787.877747" "787.652588" "787.197754" "786.970398" 
		     "786.742981" "786.515564" "786.288147" "786.060730" "785.833313" "784.696289" "784.468872" 
		     "784.241455" "784.014099" "783.786682" "783.331848" "783.104431" "782.877014" "782.649658" 
		     "782.647339" "782.422241" "782.194824" "781.967407" "781.739990" "781.512573" "781.285156" 
		     "780.830383" "780.825806" "780.750000" "780.674194" "780.602966" "780.600708" "780.598389" 
		     "780.522583" "780.446838" "780.375549" "780.371033" "780.295227" "780.219421" "780.148132" 
		     "780.145874" "780.143616" "780.029907" "779.920715" "779.916199" "779.802490" "779.693359" 
		     "779.688782" "779.575073" "779.465942" "779.461365" "779.347656" "779.238525" "779.233948" 
		     "779.120239" "779.011108" "779.006592" "778.892883" "778.783691" "778.781433" "778.779175" 
		     "778.665466" "778.556274" "778.554016" "778.551758" "778.438049" "778.328918" "778.324341" 
		     "778.210632" "778.101501" "778.099182" "778.096924" "777.983215" "777.874084" "777.871826" 
		     "777.869507" "777.785034" "777.700623" "777.646667" "777.616150" "777.531677" "777.447205" 
		     "777.362732" "777.278259" "777.193787" "777.191833" "777.109314" "777.024841" "776.964417" 
		     "776.940430" "776.855957" "776.771484" "776.737061" "776.687012" "776.602539" "776.518066" 
		     "776.509644" "776.433594" "776.349121" "776.282227" "776.275391" "776.264709" "776.180237" 
		     "776.150330" "776.095764" "776.054810" "776.025269" "776.011292" "775.926819" "775.900208" 
		     "775.842346" "775.827393" "775.775085" "775.757874" "775.673401" "775.650024" "775.599976" 
		     "775.588928" "775.524963" "775.504517" "775.420044" "775.399902" "775.372620" "775.370300" 
		     "775.335571" "775.274780" "775.251099" "775.166626" "775.149719" "775.145203" "775.082153" 
		     "775.024658" "774.997681" "774.917786" "774.913208" "774.899597" "774.774536" "774.714233" 
		     "774.690369" "774.649414" "774.524353" "774.515259" "774.462952" "774.399292" "774.316284" 
		     "774.274231" "774.235535" "774.149170" "774.117310" "774.024048" "774.008118" "773.918335" 
		     "773.898987" "773.780762" "773.773926" "773.719360" "773.649475" "773.525085" "773.520386" 
		     "773.400635" "773.325928" "773.323669" "773.321350" "773.276184" "773.261169" "773.200989" 
		     "773.151794" "773.140808" "773.080566" "773.027344" "773.020386" "772.960205" "772.902893" 
		     "772.900024" "772.871094" "772.839783" "772.779602" "772.778503" "772.719421" "772.659241" 
		     "772.654053" "772.643677" "772.598999" "772.538818" "772.529602" "772.478638" "772.418457" 
		     "772.416321" "772.405212" "772.358215" "772.298035" "772.280762" "772.237854" "772.188904" 
		     "772.177673" "772.156311" "772.117432" "772.057251" "772.031921" "771.997070" "771.961487" 
		     "771.936890" "771.907471" "771.876648" "771.816467" "771.783020" "771.756287" "771.734070" 
		     "771.696106" "771.658630" "771.635864" "771.575684" "771.534180" "771.515503" "771.506653" 
		     "771.455322" "771.409729" "771.395081" "771.334900" "771.285339" "771.279236" "771.274719" 
		     "771.160889" "771.036438" "770.911987" "770.824463" "770.787598" "770.663147" "770.597046" 
		     "770.538696" "770.414307" "770.369629" "770.289856" "770.165405" "770.142212" "770.041016" 
		     "769.916565" "769.914795" "769.792114" "769.687378" "769.667725" "769.543274" "769.460022" 
		     "769.418823" "769.294434" "769.232605" "769.169983" "769.045532" "769.005188" "768.921143" 
		     "768.796692" "768.777771" "768.672241" "768.550354" "768.547852" "768.423401" "768.322937" 
		     "768.298950" "768.174561" "768.095581" "768.050110" "767.925659" "767.868164" "767.801270" 
		     "767.676819" "767.640747" "767.552368" "767.427979" "767.413330" "767.303528" "767.185913" 
		     "767.179077" "766.958496" "766.731140" "766.503723" "766.276306" "766.048889" "765.821472" 
		     "765.594055" "765.366638" "765.364380" "765.139282" "764.911865" "764.684448" "764.682190" 
		     "764.457031" "764.229614" "764.002197")

# Calculate number of jobs to run
njobs=$((lat_max-lat_min))
njobs=$((njobs/delta+1))
njobs=$((njobs*${#time_arr[@]}))

SECONDS=0

for time in "${time_arr[@]}"; do
    # Put $time into "time_target.txt" file, which serves as input to the job run below
    magneticfield_str="T"$time 
    echo -n "" > file.log
    echo -n "$time" > time_target.txt
    echo $time
     
    outfile_prefix="/mnt/sdb/jsv/big_data/ionization_project/rig_files/"$magneticfield_str"_"$ddir"_lat_"
    for lat in $(seq $lat_min $delta $lat_max); do
        # Generate inputfile and replace lat value
        n_car=${#lat}
        if [ $n_car = 1 ]; then
            new="    "$lat
        fi
        if [ $n_car = 2 ]; then
            new="   "$lat
        fi
        if [ $n_car = 3 ]; then
            new="  "$lat
        fi
 
	# Replace lat in dummy inputfile and make new file
        search_string="%s/"$original_lat"/"$new"/g|x"
        infile=$inputfile_root"_lat_"$lat"_"$ddir 
        cp $dummy_input $infile               
        ex -s -c "$search_string" $infile     
                                              
        # Generate output filename
        outfile=$outfile_prefix$lat".txt"
        #echo -n $outfile >> file.log
        # Run model
        ./$executable $outfile $infile >> file.log & 
        let njob++
        duration=$SECONDS
	echo "Running job #("$njob"/"$njobs") for lat="$lat". Output: "$outfile". Elapsed time: $(($duration / 3600)) hours, $((($duration%3600)/60)) minutes. Estimated time left: $((($njobs-$njob)/$njob*$duration/3600)) hours, $((((($njobs-$njob)/$njob*$duration)%3600)/60)) minutes."
        # Only start new job if number of running jobs is below 30
        while [ $(pgrep -c ${executable%.*}) -gt $nmax ] ;do
            sleep 3
        done    
    done
done
echo "SCRIPT: DONE"
