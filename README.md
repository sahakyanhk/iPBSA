# iPBSA

iPBSA is a script for docking results rescoring, it minimizes docked receptor-ligand conformations in implicit solvent and calculates the binding free energy with MM/PB(GB)SA methods. The algorithm is based on a freely-available [AmberTools18](https://ambermd.org/AmberTools.php) which can be easily installed via conda  

**Usage** ./iPBSA.sh -r {.pdb} -l {dir} -n int -c bcc | gas
	
-r receptor file in pdb format  
-l path to the directory with small molecules (ligands should be in pdb format with added hydrogens)  
-n number of parallel threads (default is 8)  
-c partial charge calculation method BCC (default) or GAS  


For more details about iPBSA, please see [**Improving virtual screening results with MM/GBSA and MM/PBSA rescoring**](https://link.springer.com/article/10.1007/s10822-021-00389-3)  
