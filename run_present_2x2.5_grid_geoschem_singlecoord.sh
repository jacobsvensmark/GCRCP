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
declare lat_nospace
declare lon_nospace
declare -i nmax=30
declare -i njob=0
declare -i njobs=35

inputfile_root="TAPE1"
original_lat="AAA.AA"
original_lon="OOOO.OO"
dummy_input="input_tape_templates/TAPE1_dummy_lon_lat_vertical"
ddir="VERTICAL"
executable="tji95_present_B.out"
magneticfield="present"
indir="input_tapes/"

lon_array=("-180.00" "-177.50" "-175.00" "-172.50" "-170.00" "-167.50" "-165.00" "-162.50" 
           "-160.00" "-157.50" "-155.00" "-152.50" "-150.00" "-147.50" "-145.00" "-142.50"
           "-140.00" "-137.50" "-135.00" "-132.50" "-130.00" "-127.50" "-125.00" "-122.50"
           "-120.00" "-117.50" "-115.00" "-112.50" "-110.00" "-107.50" "-105.00" "-102.50"
           "-100.00" " -97.50" " -95.00" " -92.50" " -90.00" " -87.50" " -85.00" " -82.50"
           " -80.00" " -77.50" " -75.00" " -72.50" " -70.00" " -67.50" " -65.00" " -62.50"
           " -60.00" " -57.50" " -55.00" " -52.50" " -50.00" " -47.50" " -45.00" " -42.50"
           " -40.00" " -37.50" " -35.00" " -32.50" " -30.00" " -27.50" " -25.00" " -22.50"
           " -20.00" " -17.50" " -15.00" " -12.50" " -10.00" "  -7.50" "  -5.00" "  -2.50"
           "   0.00" "   2.50" "   5.00" "   7.50" "  10.00" "  12.50" "  15.00" "  17.50"
           "  20.00" "  22.50" "  25.00" "  27.50" "  30.00" "  32.50" "  35.00" "  37.50"
           "  40.00" "  42.50" "  45.00" "  47.50" "  50.00" "  52.50" "  55.00" "  57.50"
           "  60.00" "  62.50" "  65.00" "  67.50" "  70.00" "  72.50" "  75.00" "  77.50"
           "  80.00" "  82.50" "  85.00" "  87.50" "  90.00" "  92.50" "  95.00" "  97.50"
           " 100.00" " 102.50" " 105.00" " 107.50" " 110.00" " 112.50" " 115.00" " 117.50"
           " 120.00" " 122.50" " 125.00" " 127.50" " 130.00" " 132.50" " 135.00" " 137.50"
           " 140.00" " 142.50" " 145.00" " 147.50" " 150.00" " 152.50" " 155.00" " 157.50"
           " 160.00" " 162.50" " 165.00" " 167.50" " 170.00" " 172.50" " 175.00" " 177.50")

lat_array=("-89.50"  "-88.00"  "-86.00"  "-84.00"  "-82.00"  "-80.00"  "-78.00"  "-76.00" 
           "-74.00"  "-72.00"  "-70.00"  "-68.00"  "-66.00"  "-64.00"  "-62.00"  "-60.00"  
           "-58.00"  "-56.00"  "-54.00"  "-52.00"  "-50.00"  "-48.00"  "-46.00"  "-44.00"  
           "-42.00"  "-40.00"  "-38.00"  "-36.00"  "-34.00"  "-32.00"  "-30.00"  "-28.00" 
           "-26.00"  "-24.00"  "-22.00"  "-20.00"  "-18.00"  "-16.00"  "-14.00"  "-12.00" 
           "-10.00"  " -8.00"  " -6.00"  " -4.00"  " -2.00"  "  0.00"  "  2.00"  "  4.00" 
           "  6.00"  "  8.00"  " 10.00"  " 12.00"  " 14.00"  " 16.00"  " 18.00"  " 20.00" 
           " 22.00"  " 24.00"  " 26.00"  " 28.00"  " 30.00"  " 32.00"  " 34.00"  " 36.00" 
           " 38.00"  " 40.00"  " 42.00"  " 44.00"  " 46.00"  " 48.00"  " 50.00"  " 52.00" 
           " 54.00"  " 56.00"  " 58.00"  " 60.00"  " 62.00"  " 64.00"  " 66.00"  " 68.00" 
           " 70.00"  " 72.00"  " 74.00"  " 76.00"  " 78.00"  " 80.00"  " 82.00"  " 84.00" 
           " 86.00"  " 88.00"  " 89.50")

nlon="${#lon_array[@]}"
nlat="${#lat_array[@]}"
njobs=$((nlon * nlat))

outfile_prefix="dat/2x2.5_grid_model_output/"$magneticfield"_"$ddir"_lat_"
for lat in "${lat_array[@]}"; do
    for lon in "${lon_array[@]}"; do
        # Generate inputfile and replace lat value
        #n_car=${#lat}
        #if [ $n_car = 1 ]; then
        #    new="    "$lat
        #fi
        #if [ $n_car = 2 ]; then
        #    new="   "$lat
        #fi
        #if [ $n_car = 3 ]; then
        #    new="  "$lat
        #fi
        #if [ $n_car = 4 ]; then
        #    new=" "$lat
        #fi
        #if [ $n_car = 5 ]; then
        #    new=""$lat
        #fi
        lat_nospace="$(echo -e "${lat}" | tr -d '[:space:]')"
        lon_nospace="$(echo -e "${lon}" | tr -d '[:space:]')"
    
        #Replace lat in dummy inputfile and make new file
        search_string_lat="%s/"$original_lat"/"$lat"/g|x"
        search_string_lon="%s/"$original_lon"/"$lon"/g|x"
        infile=$indir$inputfile_root"_lat_"$lat_nospace"_lon_"$lon_nospace"_"$ddir 
        cp $dummy_input $infile               
	#echo $infile
	#echo $search_string_lat
        ex -s -c "$search_string_lat" $infile
        ex -s -c "$search_string_lon" $infile
 
        # Generate output filename
        outfile=$outfile_prefix$lat_nospace"_lon_"$lon_nospace".txt"
        logfile="logs"/$magneticfield"_"$ddir"_lat_"$lat_nospace".log"
        echo -n $outfile >> $logfile
       
        # Run model
        ./$executable $outfile $infile >> $logfile & 
        let njob++
        echo "Running job number "$njob" of "$njobs" for lat="$lat_nospace" and lon="$lon_nospace". Output in file: "$outfile
     
        # Only start new job if number of running jobs is below 30
        #while [ $(ps -C $executable | wc -l) -gt $nmax ] ;do sleep 3;done
        while [ $(pgrep -c ${executable%.*}) -gt $nmax ] ;do
            sleep 3
        done    
    done
done
echo "SCRIPT: DONE"
