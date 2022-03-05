#!/bin/bash

files=("$@")

#mkdir -m 777 ${files[0]}_MIR
#mkdir -m 777 ${files[0]}_MIRIAD
#mkdir -m 777 ${files[0]}_CASA

#/opt/rsi/idl_6.2/bin/idl -e ".r ./mir_output_to_miriad.pro" -args ${files[@]}

i=1
for element in ${@:2}
do
    asic_presence=false
    asic="asic_"
    swarm="swarm_"
    directory=${files[0]}"_MIR"/*
    pilot_data=true

    if [ "$pilot_data" == true ] : then
    then
        if [[ -n $(find ./${files[0]}_MIR -maxdepth 1 -name "*.if1.*") ]]
        then
            lsb=${files[0]}_p${i}.lsb.if1
            usb=${files[0]}_p${i}.usb.if1
            echo $lsb
            echo $usb

            uvflag vis=${files[0]}_MIR/$lsb.miriad edge=8,8,0 flagval=flag
            uvflag vis=${files[0]}_MIR/$usb.miriad edge=8,8,0 flagval=flag

            for a in $(seq 48); do
                echo $a
                fits in=${files[0]}_MIR/$lsb.miriad op=uvout select="window($a)" line=channel out=${files[0]}_MIRIAD/$lsb.$a.fits
                fits in=${files[0]}_MIR/$usb.miriad op=uvout select="window($a)" line=channel out=${files[0]}_MIRIAD/$usb.$a.fits
            done
        elif [[ -n $(find ./${files[0]}_MIR -maxdepth 1 -name "*.if2.*") ]]
        then
            lsb=${files[0]}_p${i}.lsb.if2
            usb=${files[0]}_p${i}.usb.if2
            echo $lsb
            echo $usb

            uvflag vis=${files[0]}_MIR/$lsb.miriad edge=8,8,0 flagval=flag
            uvflag vis=${files[0]}_MIR/$usb.miriad edge=8,8,0 flagval=flag

            for a in $(seq 48); do
                echo $a
                fits in=${files[0]}_MIR/$lsb.miriad op=uvout select="window($a)" line=channel out=${files[0]}_MIRIAD/$lsb.$a.fits
                fits in=${files[0]}_MIR/$usb.miriad op=uvout select="window($a)" line=channel out=${files[0]}_MIRIAD/$usb.$a.fits
            done
        else
            lsb=${files[0]}_p${i}.lsb
            usb=${files[0]}_p${i}.usb
            echo $lsb
            echo $usb

            uvflag vis=${files[0]}_MIR/$lsb.miriad edge=8,8,0 flagval=flag
            uvflag vis=${files[0]}_MIR/$usb.miriad edge=8,8,0 flagval=flag

            for a in $(seq 48); do
                echo $a
                fits in=${files[0]}_MIR/$lsb.miriad op=uvout select="window($a)" line=channel out=${files[0]}_MIRIAD/$lsb.$a.fits
                fits in=${files[0]}_MIR/$usb.miriad op=uvout select="window($a)" line=channel out=${files[0]}_MIRIAD/$usb.$a.fits
            done
        fi

    if [[ -n $(find ./${files[0]}_MIR -maxdepth 1 -name "*_${asic}${i}.*") ]]
    then
        echo $asic_presence
        lsb_asic=${files[0]}.lsb_asic_$i
        usb_asic=${files[0]}.usb_asic_$i

        uvflag vis=${files[0]}_MIR/$lsb_asic.miriad edge=8,8,0 flagval=flag
        uvflag vis=${files[0]}_MIR/$usb_asic.miriad edge=8,8,0 flagval=flag

        for a in $(seq 48); do
            echo $a
            fits in=${files[0]}_MIR/$lsb_asic.miriad op=uvout select="window($a)" line=channel out=${files[0]}_MIRIAD/$lsb_asic.$a.fits
            fits in=${files[0]}_MIR/$usb_asic.miriad op=uvout select="window($a)" line=channel out=${files[0]}_MIRIAD/$usb_asic.$a.fits
        done

        asic_presence=true
        echo $asic_presence
    fi

    if [[ -n $(find ./${files[0]}_MIR -maxdepth 1 -name "*_${swarm}${i}.*") ]]
    then

        lsb_swarm=${files[0]}.lsb_swarm_$i
        usb_swarm=${files[0]}.usb_swarm_$i

        uvflag vis=${files[0]}_MIR/$lsb_swarm.miriad edge=100,100,0 flagval=flag, log=log_file_lsb.txt
        uvflag vis=${files[0]}_MIR/$usb_swarm.miriad edge=100,100,0 flagval=flag, log=log_file_usb.txt

        if  [ "$asic_presence" = true ] ; then
            for a in `seq 49 52`; do
                echo $a
                fits in=${files[0]}_MIR/$lsb_swarm.miriad op=uvout select="window($a)" line=channel out=${files[0]}_MIRIAD/$lsb_swarm.$a.fits
                fits in=${files[0]}_MIR/$usb_swarm.miriad op=uvout select="window($a)" line=channel out=${files[0]}_MIRIAD/$usb_swarm.$a.fits
            done
        else
            for a in `seq 1 4`; do
                echo $a
                fits in=${files[0]}_MIR/$lsb_swarm.miriad op=uvout select="window($a)" out=${files[0]}_MIRIAD/$lsb_swarm.$a.fits
                fits in=${files[0]}_MIR/$usb_swarm.miriad op=uvout select="window($a)" out=${files[0]}_MIRIAD/$usb_swarm.$a.fits
            done
        fi
    fi
    i=$(($i+1))
done

#sourcename=${files[0]}

#/reduction/czdata/final/final_images/Imaging_script/anaconda2/bin/python2.7 ./phasecenter.py $sourcename
#j2000=$(<j2000.txt)
#rm j2000.txt

#/reduction/czdata/final/final_images/Imaging_script/casa-release-5.1.1-5.el6/bin/casa --nologger --nologfile -c ./CASA_tclean_continuum_only.py $j2000 ${files[@]}

#/opt/casa-release-5.1.2-4.el6/bin/mpicasa -n 3 /opt/casa-release-5.1.2-4.el6/bin/casa --nogui -c ./CASA_tclean_continuum_only.py $j2000 ${files[@]}

