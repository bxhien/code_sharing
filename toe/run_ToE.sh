# Purpose: Generate a ToE nc files with variables as follows: 
#         (1) Mean of time series of input_var within current climatology (here 1861-2005)
#         (2) Stddev of time series of input_var within current climatology (here 1861-2005)
#         (3) Upper bound as the critial value ( mean of input_var + 2*(1) )
#         (4) Time series of consecutive_year-year running mean of input_var
#         * Note that all the above variables are function of lat/lon (i.e., 2 dimensional)
#         (5) Time of Emergence (ToE) when consecutive N years exceed the upper bound (3-D, including consecutive_year-year dimension)

# Author: Yi-Xian Li, modified by Hien Bui
# Time: 2021/08/03

#---------------------------------------------------------
# Make sure the 'dir' here and the ('dir'/'diro') settings in the 2 NCL files
# are pointing to the correct directories
dir='./'

N=5 # every N year (your desired output N-year interval)

start_n=51 # starting model number
end_n=100 # ending model number
# You can run from 0 to 99 now as the limit of threads are determined by
# how large your MEM is

#---------------------------------------------------------
# Default below, basically no need to modify (at least for CESM2L) unless you change the ncl names
nclname='test.ncl'

#---------------------------------------------------------
# Processing core below

for n in $(seq $start_n 1 $end_n); # current_nmodel
do
  echo "n="$n". Thread processing years:" $i"-"$j
  nohup ncl current_nmodel=$n N=$N $dir$nclname &
# nohup stands for "ignore the HUP (hangup) signal", i.e., it'd run till finished or killed
done
wait
