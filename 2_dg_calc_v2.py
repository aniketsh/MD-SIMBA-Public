#!/home/magarkar/anaconda3/bin/python
import parmed as pmd
import subprocess, os
import MDAnalysis as mda

#1 Process gromacs trajectory
#print("Step 1/6 Postprocessing trajectory")

#subprocess.getoutput("echo 0 | gmx trjconv -s prod.tpr -f prod.xtc -o whole.xtc -pbc whole -e 1000")
#subprocess.getoutput("echo 0 | gmx trjconv -s prod.tpr -f whole.xtc -o nojump.xtc -pbc nojump")
#subprocess.getoutput("echo Protein System | gmx trjconv -s prod.tpr -f nojump.xtc -o energy.xtc -pbc mol -center")
#subprocess.getoutput("echo Protein System | gmx trjconv -s prod.tpr -f nojump.xtc -o energy.gro -e 0 -pbc mol -center")

#2 parmed to convert gromacs to AMBER
print("Step 2/6 Convering to AMBER topology")

top = pmd.load_file("../topol2.top",defines=dict(DEFINE='DUMMY'))
struct = pmd.load_file("traj.gro")

top.save('complex.parm7',overwrite=True)
struct.save('complex.rst7',overwrite=True)

#3 MDAnalysis xtc to netcdf
print("Step 3/6 Convering to AMBER trajectory")

u = mda.Universe('traj.gro', 'traj.xtc')
u.atoms.write('energy.ncdf', frames='all')

#4 Prepare trajectory with waters
print("Step 4/6 Preparing trajectory with waters")

subprocess.getoutput("cpptraj -p complex.parm7 -i /data/projects/SIMBA_PROJECT/scripts/dg_calc/closest20.in")

#5 Separate components
print("Step 5/6 Separating topology components")

subprocess.getoutput('ante-MMPBSA.py -p N20.complex.parm7 -c complex_nowat.prmtop -r protein.prmtop -l ligand.prmtop -s ":NA,CL,Na+,Cl-" -n ":MOL" --radii=mbondi2') 

#6 Calculate DG
print("Step 6/6 Running MMGBSA")

subprocess.getoutput("MMPBSA.py -O -i /data/projects/SIMBA_PROJECT/scripts/dg_calc/mmbpsa.in -o result.dat -sp N20.complex.parm7 -cp complex_nowat.prmtop -rp protein.prmtop -lp ligand.prmtop -y N20.nc")

subprocess.getoutput("rm _MM* *#*")
#subprocess.getoutput("rm whole.xtc nojump.xtc energy.ncdf N20.nc")

u2 = mda.Universe("N20.complex.parm7","N20.nc")
u2.atoms.write('for_simba.xtc', frames='all')
u2.atoms.write('for_simba.gro', frames=[0])

subprocess.getoutput("/home/magarkar/anaconda3/envs/gromacs/bin/gmx editconf -f for_simba.gro -o for_simba.pdb")
