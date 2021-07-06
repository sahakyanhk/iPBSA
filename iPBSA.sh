#!/bin/bash


display_usage() {
    echo -e "\niPBSA is a script for docking results rescoring"
    echo -e "\nUsage: ./iPBSA.sh -r <.pdb> -l <dir> -n <int> -c <bcc/gas> \n"
}

display_help() {
    echo -e "\niPBSA is a script for docking results rescoring"
    echo -e "\nUsage: ./iPBSA.sh -r <.pdb> -m <dir> -n <int> -c <bcc/gas> \n"
    echo "-r receptor file in pdb format"
    echo "-l path to the directory with ligands"
    echo "-n nukber of parallel threads (defolt is 8)" 
    echo -e "-c charge method BCC (defolt) or GAS \n"
}

if [  $# -le 1 ]
then
    display_usage
    exit 1
fi


if [[ ( $# == "--help")  ||  $# == "-h" ]]
then
    display_help
    exit 0
fi




charge=bcc
N=8

while getopts r:l:n:c: flag
do
    case "${flag}" in
        r) receptor=${OPTARG};;
        l) mols=${OPTARG};;
        n) N_thr=${OPTARG};;
        c) charge=${OPTARG};;
    esac
done

mkdir inp parm min gbsa pbsa
cp $receptor parm/rec.pdb

#set the number of cores
N=$N_thr

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
quit
EOF

cat > ./inp/prot_tleap.in <<EOF
source leaprc.protein.ff14SB
source leaprc.gaff
set default PBradii mbondi3
PROT=loadpdb parm/rec.pdb
saveamberparm PROT parm/protein.prmtop parm/protein.inpcrd
quit
/
EOF

cat > ./inp/min.in <<EOF
Initial minimisation of rec-lig complex
&cntrl
imin=1, 
maxcyc=2500, ncyc=100,
cut=16, ntb=0, igb=0,
ntpr=500,
ntwx=500,
ntwr=500
drms=0.01
ibelly=1,
bellymask=':UNL <@5'
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
igb=2,
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
for lig in ../$mols/*.pdb; do
    (
    mol=`basename $lig .pdb`
    mkdir $mol
    cd $mol
    pwd
    antechamber -i ../../$mols/${mol}.pdb -fi pdb -o $mol.mol2 -fo mol2 -at gaff -c $charge -rn UNL  -nc 0
    sleep $(( (RANDOM % 3) + 1))
    ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
    wait -n
    fi
done
wait
cd ..

#antechamber (frcmod)
cd parm
for lig in ../$mols/*.pdb
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

#tleap prot
tleap -f inp/prot_tleap.in

#tleap inputs preparation 
for lig in ./$mols/*.pdb
do
    mol=`basename $lig .pdb`
    sed "s/_mol_/$mol/g" inp/xtleap.in  > inp/tleap_$mol.in
done

#tleap complex
for lig in ./$mols/*.pdb; do
    (
    mol=`basename $lig .pdb`
    tleap -f inp/tleap_$mol.in
    sleep $(( (RANDOM % 3) + 1))
    ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
    wait -n
    fi
done

#rm inp/tleap_*.in

#minim
for lig in ./$mols/*.pdb; do
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

#GBSA
cd  gbsa
(
for lig in ../$mols/*.pdb
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


#GBSA processing
rm  rawGBSA.dat FINAL_GBSA.dat
for lig in ./$mols/*.pdb
do
    mol=`basename $lig .pdb`
    GB=`sed -e '1,/DELTA TOTAL/d' gbsa/en_pre-frames_$mol.dat | awk -F ',' '{print $NF}'|  head -n 1`
    echo $mol $GB >> rawGBSA.dat
done
sort -V rawGBSA.dat | grep ' ' > FINAL_GBSA.dat

#PBSA
cd  pbsa
(
for lig in ../$mols/*.pdb
do
    mol=`basename $lig .pdb`
    mkdir $mol
    cd $mol
    ((i=i%N)); ((i++==0)) && wait
    MMPBSA.py -O -i ../../inp/pbsa.in \
        -o  ../FINAL_RESULTS_pbsa_${mol}.dat \
        -eo ../en_pre-frames_${mol}.dat \
        -cp ../../parm/complex_${mol}.prmtop \
        -rp ../../parm/protein.prmtop \
        -lp ../../parm/${mol}.prmtop \
        -y  ../../min/min_${mol}.rst* &
    cd ..
done
)
cd ..


#PBSA processing
rm rawPBSA.dat rawPBSA.dat PBSA.dat
for lig in ./$mols/*.pdb
do
    mol=`basename $lig .pdb`
    PB=`sed -e '1,/DELTA TOTAL/d' pbsa/en_pre-frames_$mol.dat | awk -F ',' '{print $NF}' | head -n 1`
    echo $mol $PB >> rawPBSA.dat
done
sort -V rawPBSA.dat | grep ' '  > FINAL_PBSA.dat


