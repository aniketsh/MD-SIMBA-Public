#!/home/magarkar/anaconda3/envs/simba/bin/python

import MDAnalysis as mda
import subprocess,os

gro='prod.gro'
xtc='prod.xtc'

u=mda.Universe(gro,xtc)

def postprocess_gromacs(xtc):
    if os.path.isfile(xtc) and os.path.isfile('prod.tpr'):
        subprocess.getoutput("echo 0 | gmx trjconv -s prod.tpr -f %s -o whole.xtc -pbc whole"%xtc)
        subprocess.getoutput("echo 0 | gmx trjconv -s prod.tpr -f whole.xtc -o nojump.xtc -pbc nojump")
        subprocess.getoutput("echo Protein System | gmx trjconv -s prod.tpr -f nojump.xtc -o center.xtc -pbc mol -center")
        subprocess.getoutput("echo C-alpha System | gmx trjconv -s prod.tpr -f center.xtc -o traj.xtc -fit rot+trans")
        subprocess.getoutput("echo C-alpha System | gmx trjconv -s prod.tpr -f center.xtc -o traj.gro -e 0 -fit rot+trans")
    else:
        print("file not found")

###############################################################################################################################


total_frames=u.trajectory.n_frames
number_of_frames=50

if total_frames>number_of_frames:
    interval=int(total_frames/number_of_frames)
    fiter = u.trajectory[::interval]
    frames = [ts.frame for ts in fiter]
    u.atoms.write('short.xtc', frames=frames)
    postprocess_gromacs("short.xtc")
else:
    postprocess_gromacs("prod.xtc")
    
if os.path.isfile("traj.xtc") and os.path.isfile("traj.gro"):
    print("Cleaning UP")
    subprocess.getoutput("rm center.xtc whole.xtc short.xtc nojump.xtc")
