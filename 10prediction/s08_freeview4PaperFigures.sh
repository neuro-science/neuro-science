#!/bin/bash
TimeNow=$(date '+%X [%x %A %W]')
echo -e "Job starts at $TimeNow ...\n\n\n"

# # # # ###########################################################################
# # # # ### II. [Journal] Figure 
# # # # ###########################################################################
# src_path="/Users/wang/Documents/FreeSurfer/fsaverage5/surf/"
# dst_path="/Users/wang/Documents/6PD/4log/45figPaper/ccc/"
# cmd_path="/Users/wang/Documents/tmp/"
# th1=0.0024 #Figure 3, 6, 10
# th2=0.012
# theHemi=(l r)
# theSide=(outer inner)
# theDolly=1
# N1=49
# N2=4
# #The hemisphere
# for (( i=0; i<2; i++ )) ; do
#   #initial position
#   (( theAzimuth=(1-i)*180 ))
#   # (( theRoll=30-i*60 ))
#   # cmd="freeview -f ${src_path}${theHemi[i]}h.inflated -cam dolly ${theDolly} azimuth ${theAzimuth} elevation -20 roll ${theRoll}"
#   cmd="freeview -f ${src_path}${theHemi[i]}h.inflated -cam dolly ${theDolly} azimuth ${theAzimuth} elevation 0 roll 0"
#   echo $cmd > ${cmd_path}tmpcmd.txt
#   #The activation loop
#   for fname in $( ls ${src_path}*${theHemi[i]}.act ); do
#     # smooth
# 	 sname=sm${fname: ${N1}:(${#fname}-${N2}-${N1})}
# 	 echo ${fname}
# 	 echo ${sname}
#     # echo ${fname}
#     mri_surf2surf --hemi ${theHemi[i]}h \
#       --s fsaverage5 \
#       --nsmooth-in 30 \
#       --src_type curv \
#       --trg_type curv \
#       --sval ${fname} \
#       --tval ${src_path}${sname}
#     # rename
#     mv ${src_path}${theHemi[i]}h.${sname} ${src_path}${sname}
#     # show the activation and save
#     cmd="freeview -f ${src_path}${theHemi[i]}h.inflated:overlay=${src_path}${sname}:overlay_threshold=${th1},${th2}:overlay_method=linear"
#     echo $cmd >> ${cmd_path}tmpcmd.txt
#     #The side of the hemisphere
#     for (( j=0; j<2; j++ )) ; do
#       (( theAzimuth=(i!=j)*180 ))
#       # (( theElevation=20-40*j ))
#       # (( theRoll=30-(i==j)*60 ))
#       # cmd="freeview -cam azimuth 180 -ss ${dst_path}${sname}_${theHemi[i]}h_${theSide[j]}.png"
#       cmd="freeview -cam azimuth 180 -ss ${dst_path}${sname}_${theHemi[i]}h_${theSide[j]}.png"
#       echo $cmd >> ${cmd_path}tmpcmd.txt
#     done
#   done
#   #quit after all output
#   echo "freeview -quit" >> ${cmd_path}tmpcmd.txt
#   # run the command
#   freeview -cmd ${cmd_path}tmpcmd.txt
# done
# TimeNow=$(date '+%X [%x %A %W]')
# echo -e "\n\n\nJob ends at $TimeNow ...\n"

# # # # ###########################################################################
# # # # ### I. [Journal] Figure 
# # # # ###########################################################################
# src_path="/Users/wang/Documents/FreeSurfer/fsaverage5/surf/"
# dst_path="/Users/wang/Documents/6PD/4log/45figPaper/ccc/"
# cmd_path="/Users/wang/Documents/tmp/"
# th1=0 #Figure 4, 9, 11
# th2=0.012
# theHemi=(l r)
# theSide=(outer inner)
# theDolly=1
# N1=49
# N2=4
# #The hemisphere
# for (( i=0; i<2; i++ )) ; do
#   #initial position
#   (( theAzimuth=(1-i)*180 ))
#   # (( theRoll=30-i*60 ))
#   # cmd="freeview -f ${src_path}${theHemi[i]}h.inflated -cam dolly ${theDolly} azimuth ${theAzimuth} elevation -20 roll ${theRoll}"
#   cmd="freeview -f ${src_path}${theHemi[i]}h.inflated -cam dolly ${theDolly} azimuth ${theAzimuth} elevation 0 roll 0"
#   echo $cmd > ${cmd_path}tmpcmd.txt
#   #The activation loop
#   for fname in $( ls ${src_path}*${theHemi[i]}.act ); do
#     # smooth
# 	 sname=sm${fname: ${N1}:(${#fname}-${N2}-${N1})}
# 	 echo ${fname}
# 	 echo ${sname}
#     # echo ${fname}
#     mri_surf2surf --hemi ${theHemi[i]}h \
#       --s fsaverage5 \
#       --nsmooth-in 30 \
#       --src_type curv \
#       --trg_type curv \
#       --sval ${fname} \
#       --tval ${src_path}${sname}
#     # rename
#     mv ${src_path}${theHemi[i]}h.${sname} ${src_path}${sname}
#     # show the activation and save
#     cmd="freeview -f ${src_path}${theHemi[i]}h.inflated:overlay=${src_path}${sname}:overlay_threshold=${th1},${th2}:overlay_method=linear"
#     echo $cmd >> ${cmd_path}tmpcmd.txt
#     #The side of the hemisphere
#     for (( j=0; j<2; j++ )) ; do
#       (( theAzimuth=(i!=j)*180 ))
#       # (( theElevation=20-40*j ))
#       # (( theRoll=30-(i==j)*60 ))
#       # cmd="freeview -cam azimuth 180 -ss ${dst_path}${sname}_${theHemi[i]}h_${theSide[j]}.png"
#       cmd="freeview -cam azimuth 180 -ss ${dst_path}${sname}_${theHemi[i]}h_${theSide[j]}.png"
#       echo $cmd >> ${cmd_path}tmpcmd.txt
#     done
#   done
#   #quit after all output
#   echo "freeview -quit" >> ${cmd_path}tmpcmd.txt
#   # run the command
#   freeview -cmd ${cmd_path}tmpcmd.txt
# done
# TimeNow=$(date '+%X [%x %A %W]')
# echo -e "\n\n\nJob ends at $TimeNow ...\n"

