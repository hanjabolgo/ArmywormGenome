###########################################
# Plot the result of BUSCO
###########################################
base=~/Project/Armyworms/

cd $base/01_Assembly/
mkdir -p Asembly_BUSCO/summary
cd Asembly_BUSCO/
find ./ -name "short_summary.specific.*.txt" | xargs -i cp {} summary/
python ./generate_plot.py --working_directory ./summary
