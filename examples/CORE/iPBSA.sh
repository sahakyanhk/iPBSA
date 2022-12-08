#!/bin/bash
#iPBSA minimizes docked receptor-ligand conformations in implicit solvent and calculates the binding free energy with MM/PB(GB)SA methods. 
#The algorithm is based on a freely-available [AmberTools](https://ambermd.org/AmberTools.php) which can be easily installed via conda  
#For more details about iPBSA, please see https://link.springer.com/article/10.1007/s10822-021-00389-3

#set -e

display_usage() {
    echo -e "\niPBSA is a script for docking results rescoring"
    echo -e "\nUsage: ./iPBSA.sh -r <.pdb> -l <lig dir> -n <int> -c <bcc/gas> -o <output dir> \n"
    echo "-r receptor file in pdb format"
    echo "-l path to the directory with small molecules (ligands should be in pdb format with added hydrogens)."
    echo "-n number of parallel threads (The number of threads must not excead the number of ligands)." 
    echo "-o output directory" 
    echo -e "-c charge method BCC (defolt) or GAS \n"
}

if [  $# -le 1 ]
then
    display_usage
    exit 1
fi


charge=bcc
N=4

while getopts r:l:o:n:c: flag
do
    case "${flag}" in
        r) receptor=${OPTARG};;
        l) mols=${OPTARG};;
        o) output=${OPTARG};;
        n) N=${OPTARG};;
        c) charge=${OPTARG};;
    esac
done

mkdir $output; cd $output


mkdir -p inp parm min gbsa pbsa
pdb4amber -i ../$receptor -o parm/rec.pdb --nohyd --dry

#prepare input files
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
echo "running antechamber..."
cd parm
for lig in ../../$mols/*.pdb; do
    (
    mol=`basename $lig .pdb`
    mkdir -p $mol
    cd $mol
    pwd
    antechamber -i ../../../$mols/${mol}.pdb -fi pdb -o $mol.mol2 -fo mol2 -at gaff -c $charge -rn UNL  -nc 0 
    sleep $(( (RANDOM % 3) + 1))
    ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
    wait -n
    fi
done > ../antechamber.log
wait
cd ..

#antechamber (frcmod)
cd parm
for lig in ../../$mols/*.pdb
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
done> ../parmchk2.log
cd ..

#tleap prot
echo "preparing receptor..."
tleap -f inp/prot_tleap.in > prot_tleap.log 

#tleap mols 
echo "preparing small molecules..."
for lig in ./../$mols/*.pdb
do
    mol=`basename $lig .pdb`
    sed "s/_mol_/$mol/g" inp/xtleap.in  >> inp/tleap_$mol.in
done

#tleap complex
echo "preparing complexes..."

for lig in ./../$mols/*.pdb; do
    (
    mol=`basename $lig .pdb`
    tleap -f inp/tleap_$mol.in 
    sleep $(( (RANDOM % 3) + 1))
    ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
    wait -n
    fi
done > tleap.log

#rm inp/tleap_*.in

#minim
echo "minimization..."
for lig in ./../$mols/*.pdb; do
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
done > minimization.log
wait

#GBSA
echo "running mmgbsa..."

cd  gbsa
(
for lig in ../../$mols/*.pdb
do
    mol=`basename $lig .pdb`
    mkdir -p $mol
    cd $mol
    ((i=i%N)); ((i++==0)) && wait
    MMPBSA.py -O -i ../../inp/gbsa.in \
        -o  ../FINAL_RESULTS_GBSA_${mol}.dat \
        -eo ../en_pre-frames_${mol}.dat \
        -cp ../../parm/complex_${mol}.prmtop \
        -rp ../../parm/protein.prmtop \
        -lp ../../parm/${mol}.prmtop \
        -y  ../../min/min_${mol}.rst*  &
    cd ..
done 

date
wait
date
) > ../gbsa.log


cd ..


#PBSA
echo "running mmpbsa..."

cd  pbsa
(

for lig in ../../$mols/*.pdb
do

    mol=`basename $lig .pdb`
    mkdir -p $mol
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

date
wait
date
) > ../pbsa.log



cd ..



for lig in ./../$mols/*.pdb
do
   mol=`basename $lig .pdb`
   GB=`sed -e '1,/DELTA TOTAL/d' gbsa/en_pre-frames_${mol}.dat | awk -F ',' '{print $NF}' |  head -n 1`
   PB=`sed -e '1,/DELTA TOTAL/d' pbsa/en_pre-frames_${mol}.dat | awk -F ',' '{print $NF}' | head -n 1`
   echo $mol $PB >> rawPBSA.dat
   echo $mol $GB >> rawGBSA.dat
done
sort -V rawGBSA.dat | grep ' ' > FINAL_GBSA.dat
sort -V rawPBSA.dat | grep ' ' > FINAL_PBSA.dat
