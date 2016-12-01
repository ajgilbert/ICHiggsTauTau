#!/bin/bash

if [ "$#" -ne "2" ]; then
    echo "Usage: $0 <doSubmit> <do4params>"
    exit 0
fi


DOSUBMIT=$1
DO4PARAMS=$2
infolder=output_run2ana_161121_ICHEP_NLO_cut/
outfolder=cards_run2ana_161121_ICHEP_NLO_cut/
blind=true
#zvvstat=18
mkdir -p $outfolder

extraoptions="--do_ues=false" #--do_ggh=false --do_separate_qcdewk=false"

for channel in qcd nunu enu munu mumu taunu #qcd #topl topb
do
    echo " ********************************"
    echo " *** Processing channel $channel"
    echo " ********************************"

    mkdir -p $outfolder/$channel

    #for mindphicut in 2.31 #1.01 1.21 1.41 1.61 1.81 2.01 2.11 2.21 2.31 2.41 2.51 2.61 2.81
    #do
    mindphicut=2.3
    if [ "$channel" == "qcd" ]; then
	mindphicut=-1
    elif [ "$channel" == "taunu" ]; then
	mindphicut=-1
    else
	mindphicut=2.3
    fi
    echo "channel $channel, mindphicut: $mindphicut"
	
    for minjjcut in 1101 #1601 #801 901 1001 1101 1201 1301 1401 1501 1601 1701 1801 1901
    do
	OUTNAME=$outfolder/$channel/vbfhinv_${channel}_13TeV_${mindphicut}_${minjjcut}.txt
	if (( "$DO4PARAMS" == "1" )); then
	    OUTNAME=$outfolder/$channel/vbfhinv_${channel}_13TeV_${mindphicut}_${minjjcut}_4params.txt
	fi
	if (( "$DOSUBMIT" == "0" )); then
	    echo "./bin/makeCountingCard -i $infolder --blind=$blind -o $OUTNAME -m 125 --channel $channel --do_latex true --do_datatop false --zvvstat 0 --qcdrate 0 --mcBkgOnly=true --do_run2=true --do_4params=$DO4PARAMS --minvarXcut=$minjjcut --minvarYcut=$mindphicut --histoToIntegrate=alljetsmetnomu_mindphi:dijet_M $extraoptions | tee $outfolder/$channel/card_${mindphicut}_${minjjcut}.log"
	else
	    ./bin/makeCountingCard -i $infolder --blind=$blind -o $OUTNAME -m 125 --channel $channel --do_latex true --do_datatop false --zvvstat 0 --qcdrate 0 --mcBkgOnly=true --do_run2=true --do_4params=$DO4PARAMS --minvarXcut=$minjjcut --minvarYcut=$mindphicut --histoToIntegrate=alljetsmetnomu_mindphi:dijet_M $extraoptions | tee $outfolder/$channel/card_${mindphicut}_${minjjcut}.log
	fi
    done
#done
done

