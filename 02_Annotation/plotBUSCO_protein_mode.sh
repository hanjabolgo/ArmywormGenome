###########################################
# Plot the result of BUSCO
###########################################
base=~/Project/Armyworms/

cd $base/02_Assembly
mkdir -p Assembly_BUSCO/summary
cd Assembly_BUSCO
find ./ -name "short_summary.specific.*.txt" | xargs -i cp {} summary/
python ./generate_plot.py --working_directory ./summary
