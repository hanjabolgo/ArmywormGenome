################################################################################
# Genome Assembly
################################################################################

################################################################################
# 1 Genome survey
################################################################################

base=~/Project/Armyworms/01_Assembly/GenomeSurvey
rawdata=~/Project/Armyworms/00_RawData/GenomeSurvey
software=~/biosoft
scripts=~/scripts/

cd $base/
mkdir MsepF MsepM MlorF MlorM

################################################################################
# 1.1 Using GCE to estimating the genome size of female M. separata
################################################################################
cd $base/MsepF
find $rawdata/MsepF | xargs ls -d > MsepF.lib
$software/gce-1.0.2/kmerfreq -k 19 -t 16 -p MsepF -w 1 -c 1 MsepF.lib
KmerIndNum=`grep -oP "(?<=#Kmer indivdual number: ).*" MsepF.kmer.freq.stat`
less MsepF.kmer.freq.stat | perl -ne 'next if(/^#/ || /^\s/); print; ' | awk '{print $1"\t"$2}' > MsepF.kmer.freq.stat.2colum
$software/gce-1.0.2/gce -f MsepF.kmer.freq.stat.2colum -g $KmerIndNum -H 1 -c 115 > MsepF_gce2.table 2> MsepF_gce2.log

################################################################################
# 1.2 Using GCE to stimating the genome size of male M. separata
################################################################################
cd $base/MsepM
find $rawdata/MsepM | xargs ls -d > MsepM.lib
$software/gce-1.0.2/kmerfreq -k 19 -t 16 -p MsepM -w 1 -c 1 MsepM.lib
KmerIndNum=`grep -oP "(?<=#Kmer indivdual number: ).*" MsepM.kmer.freq.stat`
less MsepM.kmer.freq.stat | perl -ne 'next if(/^#/ || /^\s/); print; ' | awk '{print $1"\t"$2}' > MsepM.kmer.freq.stat.2colum
$software/gce-1.0.2/gce -f MsepF.kmer.freq.stat.2colum -g $KmerIndNum -H 1 -c 68 > MsepF_gce2.table 2> MsepF_gce2.log

################################################################################
# 1.3 Using GCE to stimating the genome size of female M. loreyi
################################################################################
cd $base/MlorF
find $rawdata/MlorF | xargs ls -d > MlorF.lib
$software/gce-1.0.2/kmerfreq -k 19 -t 16 -p MlorF -w 1 -c 1 MlorF.lib
KmerIndNum=`grep -oP "(?<=#Kmer indivdual number: ).*" MlorF.kmer.freq.stat`
less MlorF.kmer.freq.stat | perl -ne 'next if(/^#/ || /^\s/); print; ' | awk '{print $1"\t"$2}' > MlorF.kmer.freq.stat.2colum
$software/gce-1.0.2/gce -f MlorF.kmer.freq.stat.2colum -g $KmerIndNum -H 1 -c 68 > MlorF_gce2.table 2> MlorF_gce2.log

################################################################################
# 1.4 Using GCE to stimating the genome size of male M. loreyi
################################################################################
cd $base/MlorM
find $rawdata/MlorM | xargs ls -d > MlorM.lib
$software/gce-1.0.2/kmerfreq -k 19 -t 16 -p MlorM -w 1 -c 1 MlorM.lib
KmerIndNum=`grep -oP "(?<=#Kmer indivdual number: ).*" MlorM.kmer.freq.stat`
less MlorM.kmer.freq.stat | perl -ne 'next if(/^#/ || /^\s/); print; ' | awk '{print $1"\t"$2}' > MlorM.kmer.freq.stat.2colum
$software/gce-1.0.2/gce -f MlorM.kmer.freq.stat.2colum -g $KmerIndNum -H 1 -c 75 > MlorM_gce2.table 2> MlorM_gce2.log
