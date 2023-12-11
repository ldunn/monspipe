#!/usr/bin/env bash

echo $(which pat)

xprof='/fred/oz002/ldunn/UTMOST2D/code/profile_nt'
example='/fred/oz002/ldunn/UTMOST2D/code/example_nt'
make_plots='/fred/oz002/ldunn/meertime_dataportal/pipeline/make_plots.sh'

TIMS='/fred/oz002/ldunn/meertime_dataportal/data/tims'
PARS='/fred/oz002/ldunn/meertime_dataportal/data/pars'
FIX_PARS='/fred/oz002/ldunn/meertime_dataportal/data/fix_pars'
PROFILES='/fred/oz002/ldunn/meertime_dataportal/data/profiles'
PLOTS='/fred/oz002/ldunn/meertime_dataportal/data/images'
DATA='/fred/oz002/ldunn/meertime_dataportal/data/post'

psrj=$1
date=$2
echo "Doing PSR ${psrj}, observation at ${date}"
if ! test -f "${PROFILES}/${psrj}.raw.std.b"; then
    echo "Do not have profile! Exiting."
    exit 1
fi

obs=$DATA/${psrj}/${date}/${psrj}_${date}.ar
if ! test -f "$obs"; then
    echo "Could not find raw archive at ${obs}! Exiting."
    exit 1
fi

jobfs_dir=$JOBFS/$psrj/$date
mkdir -p $jobfs_dir
cp $obs $jobfs_dir

pam -e decim -pb2 $jobfs_dir/${psrj}_${date}.ar
psredit -m -c site=MOST_NS $jobfs_dir/${psrj}_${date}.decim

if test -f "$FIX_PARS/${psrj}.par"; then
    echo "Updating ephemeris for ${psrj}"
    pam -m -E $FIX_PARS/${psrj}.par $jobfs_dir/${psrj}_${date}.decim
fi

paz -m -E 5 $jobfs_dir/${psrj}_${date}.decim
paz -m -F "805 808" -F "855 858" $jobfs_dir/${psrj}_${date}.decim

$example $jobfs_dir/${psrj}_${date}.decim
$xprof -ps ${PROFILES}/${psrj}.raw.std.b $jobfs_dir/${psrj}_${date}.decim.raw -Z -z > $jobfs_dir/${psrj}_${date}_zap.psh
chmod +x $jobfs_dir/${psrj}_${date}_zap.psh
sed -i '/^DEBUG/d' $jobfs_dir/${psrj}_${date}_zap.psh
$jobfs_dir/${psrj}_${date}_zap.psh $jobfs_dir/${psrj}_${date}.decim -e clean
pam -DFT -e FT $jobfs_dir/${psrj}_${date}.clean
snr=`psrstat -c snr=standard -c snr:file=${PROFILES}/${psrj}.std -c snr $jobfs_dir/${psrj}_${date}.FT | grep -oP 'snr=\K(\d+(.\d+)?)'`
psrstat -c snr=standard -c snr:file=${PROFILES}/${psrj}.std -c snr $jobfs_dir/${psrj}_${date}.FT
pushd $jobfs_dir
pat -s ${PROFILES}/${psrj}.std -A FDM -f "tempo2 IPTA" -t -K ${psrj}_${date}_pat.png/PNG -X "-snr ${snr}" ${psrj}_${date}.FT
tim=$(pat -s ${PROFILES}/${psrj}.std -A FDM -f "tempo2 IPTA" -X "-snr ${snr}" ${psrj}_${date}.FT | grep mo)
$make_plots $psrj $date
popd
echo $tim > $DATA/$psrj/$date/toa.tim

mkdir -p $DATA/$psrj/$date
mkdir -p $PLOTS/$psrj/$date
mv $jobfs_dir/${psrj}_${date}_zap.psh $DATA/$psrj/$date/zap.psh
mv $jobfs_dir/${psrj}_${date}.clean $DATA/$psrj/$date/
mv $jobfs_dir/${psrj}_${date}.FT $DATA/$psrj/$date/
mv $jobfs_dir/${psrj}_${date}_*.png $PLOTS/$psrj/$date/
rm -r $jobfs_dir
