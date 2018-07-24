#!/bin/bash
# Use only for TSoft files created directly by iGrav program (must have the same folder structure, i.e. some_folder/iGravNNN_YYYY/MMDD/)
# Check if all iGrav files in given year (part of input PATH) contain required string (e.g. last time stamp "23 59 59")
# Will print all files without such string
# Usage: bash check_igrav_contains.sh /path/to/iGrav033_Data/iGrav033_2017 "23 59 59"

curfolder=$(pwd) # get current folder to switch back after checking data
cd "$1" # switch to folder with data
folders=$(ls -d */) # list all folders (e.g. "1028 1029"...)
# run for loop for all iGrav subFolders and all tsoft files inside
for d in $folders
do 
   cd "$d"
   for f in $(ls *tsf)
   do 
     out=$(grep "$2" "$f");
	if [ -z "$out" ]; then 
	   echo $f
	fi
   done
   cd ..
done
cd $curfolder # switch back to starting folder
