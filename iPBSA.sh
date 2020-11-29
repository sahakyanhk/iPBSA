#!/bin/bash

mkdir inp parm min gbsa pbsa
cp rec.pdb parm

#set the number of cores
N=8

cat > ./inp/xtleap.in <<EOF
source leaprc.protein.ff14SB
source leaprc.gaff
set default PBradii mbondi3
PROT=loadpdb parm/rec.pdb
LIG=loadmol2 parm/_mol_/_mol_.mol2
loadamberparams parm/_mol_/_mol_.frcmod
complex = combine {PROT LIG}
saveamberparm complex parm/complex__mol_.prmtop parm/complex__mol_.inpcrd
saveamberparm LIG parm/_mol_.prmtop parm/_mol_.inpcrd
savepdb complex parm/complex__mol_.pdb
quit
EOF

cat > ./inp/prot_tleap.in <<EOF
source leaprc.protein.ff14SB
source leaprc.gaff
set default PBradii mbondi3
PROT=loadpdb parm/rec.pdb
saveamberparm PROT parm/protein.prmtop parm/protein.inpcrd
quit
EOF

cat > ./inp/min.in <<EOF
Initial minimisation of rec-lig complex
&cntrl
imin=1, 
maxcyc=1500, ncyc=500,
cut=16, ntb=0, igb=8,
ntpr=100,
ntwx=100,
ntwr=-100
drms=0.01
&end
EOF

cat > ./inp/gpu_min.in <<EOF
Initial minimisation of rec-lig complex
&cntrl
imin=1, 
maxcyc=5000, ncyc=500,
cut=9999, ntb=0, igb=8,
ntpr=100,
ntwx=100,
ntwr=100,
drms=0.01
&end
/
EOF

cat > ./inp/gbsa.in <<EOF
mmgbsa  analysis
&general
startframe=1, interval=1, endframe=9999,
verbose=2, keep_files=0, netcdf=1,
strip_mask = ":WAT,Cl-,Na+,K+",
/
&gb
igb=8
saltcon=0.150,
/
EOF

cat > ./inp/pbsa.in <<EOF
mmgbsa  analysis
&general
startframe=1, interval=1, endframe=9999,
verbose=2, keep_files=0, netcdf=1,
strip_mask = ":WAT,Cl-,Na+,K+",
/
&pb
inp=1
istrng=0.150
/
EOF



#antechamber
cd parm
for lig in ../mols/*.pdb; do
    (
    mol=`basename $lig .pdb`
    mkdir $mol
    cd $mol
    antechamber -i ../../mols/${mol}.pdb -fi pdb -o $mol.mol2 -fo mol2 -at gaff -c bcc -rn UNL  -nc 0
    sleep $(( (RANDOM % 3) + 1))
    ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
    wait -n
    fi
done
wait
cd ..



#antechamber
cd parm
for lig in ../mols/*.pdb
do
    (
    mol=`basename $lig .pdb`
    cd $mol
    parmchk2 -i ${mol}.mol2 -f mol2 -o ${mol}.frcmod
    cd ..
    sleep $(( (RANDOM % 3) + 1))
    ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
    wait -n
    fi
done
cd ..

#   #tleap prot
tleap -f inp/prot_tleap.in

##tleap complex 
for lig in ./mols/*.pdb
do
    mol=`basename $lig .pdb`
    sed "s/_mol_/$mol/g" inp/xtleap.in  > inp/tleap_${mol}.in
done



for lig in ./mols/*.pdb; do
    (
    mol=`basename $lig .pdb`
    tleap -f inp/tleap_${mol}.in
    sleep $(( (RANDOM % 3) + 1))
    ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
    wait -n
    fi
done
#   
rm inp/tleap_*.in

#minim
for lig in ./mols/*.pdb; do
    (
    mol=`basename $lig .pdb`
    sander -O -i inp/min.in \
           -p parm/complex_${mol}.prmtop \
           -c parm/complex_${mol}.inpcrd \
           -r min/min_${mol}.rst \
           -ref parm/complex_${mol}.inpcrd \
           -o min/minim_${mol}.out 
    sleep $(( (RANDOM % 3) + 1))
    ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
    wait -n
    fi
done
wait


#   for lig in ./mols/*.pdb; do
#       mol=`basename $lig .pdb`
#       pmemd.cuda  -O -i inp/gpu_min.in \
#              -p parm/complex_${mol}.prmtop \
#              -c parm/complex_${mol}.inpcrd \
#              -ref parm/complex_${mol}.inpcrd \
#              -r min/min_${mol}.rst \
#              -o min/minim_${mol}.out 
#   done
#   



#GBSA
cd  gbsa
(
for lig in ../mols/*.pdb
do
    mol=`basename $lig .pdb`
    mkdir $mol
    cd $mol
    ((i=i%N)); ((i++==0)) && wait
    MMPBSA.py -O -i ../../inp/gbsa.in \
        -o  ../FINAL_RESULTS_GBSA_${mol}.dat \
        -eo ../en_pre-frames_${mol}.dat \
        -cp ../../parm/complex_${mol}.prmtop \
        -rp ../../parm/protein.prmtop \
        -lp ../../parm/${mol}.prmtop \
        -y  ../../min/min_${mol}.rst* &
    cd ..
done
)
cd ..


#process PBSA/GBSA results
rm GBSA.dat FINAL_GBSA.dat
for lig in ./mols/*.pdb
do
    mol=`basename $lig .pdb`
    GB=`awk '/DELTA TOTAL/,EOF' gbsa/en_pre-frames_${mol}.dat | awk -F ',' '{print $NF}'| sort -n | head -n 1`
    echo $mol,$GB >> rawGBSA.csv
done
sort -V rawGBSA.csv > GBSA.csv
grep ',[-][0-9]' GBSA.csv > FINAL_GBSA.csv

#pbsa
cd  pbsa
(
for lig in ../mols/*.pdb
do
    mol=`basename $lig .pdb`
    mkdir $mol
    cd $mol
    ((i=i%N)); ((i++==0)) && wait
    MMPBSA.py -O -i ../../inp/pbsa.in \
        -o  ../FINAL_RESULTS_pbsa_${mol}.dat \
        -eo en_for_frames.dat \
        -cp ../../parm/complex_${mol}.prmtop \
        -rp ../../parm/protein.prmtop \
        -lp ../../parm/${mol}.prmtop \
        -y  ../../min/min_${mol}.rst* &
    cd ..
done
)
cd ..

#process PBSA/pbsa results
rm pbsa.dat FINAL_pbsa.dat
for lig in ./mols/*.pdb
do
    mol=`basename $lig .pdb`
    GB=`awk '/DELTA TOTAL/,EOF' gbsa/en_pre-frames_${mol}.dat | awk -F ',' '{print $NF}'| sort -n | head -n 1`
    echo $mol,$GB >> rawPBSA.csv
done
sort -V rawPBSA.csv > PBSA.csv
grep ',[-][0-9]' PBSA.csv > FINAL_PBSA.csv


