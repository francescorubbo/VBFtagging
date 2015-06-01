#!/bin/bash

[ "$USER" == "rubbo" ] && WorkDir=/u/at/rubbo/nfs/VBFtagging/
# add similar line if you are not kurinsky

SubFileLoc=`pwd`/_submitSingleJob.sh
DateSuffix=`date +%Y%m%d_%Hh%Mmin`

export BSUB_QUIET=

echo '#!/bin/bash
echo CD to $1
echo CMD is $4

cd $1
source setup.sh
cmd=$4

echo MAKING TEMP DIR $2
JOBFILEDIR=$2
mkdir $JOBFILEDIR
REALOUT=$3
echo MADE TEMP DIR $JOBFILEDIR
echo WILL COPY TO $REALOUT

shift 4
echo Calling $cmd $*
$cmd $*
cp -r $JOBFILEDIR/*.root $REALOUT
echo COPYING to $REALOUT
rm -rf $JOBFILEDIR
' > $SubFileLoc
chmod u+x $SubFileLoc

#----------------
# OPTIONS

Process=4
mus="25"
Queue=short
nevents=20
njobs=500

OutDirFinal=`pwd`/files/${DateSuffix}
mkdir -p $OutDirFinal
echo
echo "Submitting $njobs jobs each with $nevents events to $Queue"

echo "Job Parameter Grid:"
echo "Mu  - "$mus

for mu in $mus ; do
    LogPrefix=`pwd`/logs/${DateSuffix}/${DateSuffix}_${mu}_${nevents}
    mkdir -p `dirname $LogPrefix`
    echo $LogPrefix
		    
    for (( ii=1; ii<=$njobs; ii++ )) ;  do
	OutDir=/tmp/${DateSuffix}_${ii}/
	
	bsub -q ${Queue} -R rhel60 -o $LogPrefix${ii}.log \
	    $SubFileLoc ${WorkDir} ${OutDir} ${OutDirFinal} \
	    ./VBFTagging \
	    --Pileup $mu                 \
	    --OutFile ${OutDir}/Sample_mu-${mu}_nevents-${nevents}_job-${ii}.root \
	    --Proc ${Process} \
	    --NEvents ${nevents} \
	    --Seed ${ii} \
	    done
    done
done
