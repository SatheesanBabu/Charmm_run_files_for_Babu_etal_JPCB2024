# Charmm_run_files_for_Babu_etal_JPCB2024
Files required to reproduce results reported in paper "Solution Ionic Strength Can Modulate Functional Loop Conformations in E. Coli Dihydrofolate Reductase" By Babu, C. Satheesan, 
This database contains eqilibrated .CRD files for CHARMM  simulations at ionic strengths of 0.075M, 0.15M, 0.3M and 0.75M. Also given are .RST restart files corresponding to .CRD files. 
The file names are arranged in the order PDBid_OCCLUDED/OPEN/CLOSED/_IonicStrength_EQ.CRD/.RST
For example file 6cw7_OCC_0.75M_EQ.CRD file corresponds to a simulation of DHFR structure 6cw7.pdb which is in OCCLUDED M20loop having a CaCl2 solution ionic strengh of 0.75M.  _EQ stands for equilibrated configuration (.RST stands for CHARMM restart file)
Files ending with .inp extension are CHARMM input files.  The files have form PDBid_OCC/OPN/CLO_Ionic Strengh_dynP.inp
File dhfr.seq is the residue sequence file
A fortran program for computing distance dependent electrostatic around the nose of M20 loop (see paper) is included with extension .f
dcd_loop.inp is the data file for this program
trajectory coorinate file is stored in file 'out_file' as inpput file (in uncompressed form)
Please contact Dr. C. S. Babu for any questions  (email:  satheesanbabu@gmail.com)
