# iPBSA

iPBSA is a script for docking results rescoring, it minimizes docked receptor-ligand conformations in implicit solvent and calculates the binding free energy with MM/PB(GB)SA methods. The algorithm is based on a freely-available [AmberTools](https://ambermd.org/AmberTools.php) which can be easily installed via conda  

```console
#Download and export iPBSA
wget https://github.com/sahakyanhk/iPBSA/archive/refs/heads/main.zip; unzip main.zip; chmod +x iPBSA-main/iPBSA.sh
export PATH=$(pwd)/iPBSA-main/:$PATH
iPBSA.sh
```

**Usage** ./iPBSA.sh -r <.pdb> -l <lig_dir>  -n <int\> -c <bcc/gas> -o <out_dir>  
	
-r receptor file in pdb format  
-l path to the directory with small molecules (ligands should be in pdb format with added hydrogens)  
-n number of parallel threads 
-c charge method BCC (defolt) or GAS  
-o output directory  

For more details, please see [**Improving virtual screening results with MM/GBSA and MM/PBSA rescoring**](https://link.springer.com/article/10.1007/s10822-021-00389-3)  
