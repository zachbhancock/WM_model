#!/bin/bash

# pick logfiles dir to search for errors in

#logdirtosearch=logfiles_gendivextrap_WM
#logdirtosearch=logfiles_marmapdists
#logdirtosearch=logfiles_ecoregions
#logdirtosearch=logfiles_popgenstats
logdirtosearch=logfiles_wmModel


# create output file:

logfile=log_summary.txt

if [ -f ./$logfile ]; then
    echo "$logfile already exists, removing and creating new one"
    rm ./$logfile
    touch ./$logfile
else 
    touch ./$logfile
fi




# search for errors in log files

echo 'Searching for non-case-sensitive: "error", "cannot stat", "cannot create", "quota exceeded", "unable to open", "not positive definite" and "failed" in .err and .out log files.' >> ./$logfile

if [ $(grep -riE 'error|cannot stat|cannot create|quota exceeded|unable to open|not positive definite|failed' ./$logdirtosearch | wc -l) -eq 0 ]
then
    echo -e 'No errors found in log file. \n' >> ./$logfile
else
    echo -e 'ERROR: Errors found in log files. The errors are: \n' >> ./$logfile
    grep -riE 'error|cannot stat|cannot create|quota exceeded|unable to open|not positive definite|failed' ./$logdirtosearch >> ./$logfile
fi
