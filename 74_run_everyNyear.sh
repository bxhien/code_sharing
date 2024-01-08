# Purpose: (1) Generate N-year nc files output (2-D lat/lon) under each model's output folder (will be auto-created)
#          (2) Combine the N-year nc files into a 3-D (time-series 2-D map) collection nc file.
#---------------------------------------------------------
# Make sure the 'dir' here and the ('dir'/'diro') settings in the 2 NCL files
# are pointing to the correct directories
dir='./'

N=5 # every N year (your desired output N-year interval)

start_n=81 # starting model number
end_n=100 # ending model number
# Suggesting "n" numbers at a time: (Number of your available CPU/thread processors) / ( (end_year-start_year+1) / N )
# E.g., 48 CPU processors & N=10 -> 2 models at a time (e.g., start_n=0 and end_n=1)
# E.g., 96 CPU processors & N=10 -> 4 models at a time (e.g., start_n=2 and end_n=5)
# E.g., 96 CPU processors & N=20 -> 8 models at a time (e.g., start_n=0 and end_n=7)
# Note that (HIST+SSP) of the same model can be only counted as one n instead of two,
# i.e., 100 HIST+ 100SSP -> min(start_n)=0 and max(end_n)=99 (total 100 models)
#---------------------------------------------------------
# Default below, basically no need to modify (at least for CESM2L) unless you change the ncl names
nclname='cal_std_everyNyear_main_prec.ncl'
combinencl='combine_everyNyear_main_prec.ncl'

start_year=1861 # start from this year
end_year=2100 # end until this year

years_per_nc=10 # how many years in a single input nc file (should be 10 for CESM-2L)

#---------------------------------------------------------
# Processing core below

end_i=$(($end_year-$N+1))
range_year=$(($end_year-$start_year+1))
if [ $(( $range_year % $N )) -ne 0 ]; then
  echo 'Check your end_year and start_year to make sure your N is divisible.'
  exit
fi

for n in $(seq $start_n 1 $end_n); # current_nmodel
do
  #for i in $(seq $start_year $N $end_i);
  for i in $(seq $start_year $N $end_year);
  do
    j=$(($i+$N-1))
    echo "Thread processing years:" $i"-"$j
    nohup ncl current_nmodel=$n N=$N start_year=$i end_year=$j $dir$nclname &
    # nohup stands for "ignore the HUP (hangup) signal", i.e., it'd run till finished or killed
  done
done
wait

for n in $(seq $start_n 1 $end_n);
do
  nohup ncl current_nmodel=$n N=$N $dir$combinencl &
done
wait
