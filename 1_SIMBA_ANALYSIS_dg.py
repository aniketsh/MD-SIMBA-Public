#!/home/magarkar/anaconda3/bin/python
import parmed as pmd
import subprocess, os, sys
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from numpy.linalg import norm
from rdkit import Chem
import numpy as np
import pandas as pd
from numpy import genfromtxt

def load_sdf_with_N(sdf_file, **kargs):
    ''' Loads a molecule correctly when there is a N(+)R4 in it by setting the charge of N to +1 '''
    mol_list = []
    try:
        tmp_mol = Chem.SDMolSupplier(sdf_file, removeHs=False, sanitize=False)
    except:
        eprint('Could not load mol from', sdf_file)
        return []
    for this_mol in tmp_mol:
        if not this_mol:
            eprint('Error loading Molecule')
            continue
        for atom in this_mol.GetAtoms():
            if atom.GetSymbol() == 'N':
                # print(atom.GetSymbol(), atom.GetExplicitValence(), atom.GetFormalCharge())
                if atom.GetExplicitValence() == 4 and atom.GetFormalCharge() == 0:
                    atom.SetFormalCharge(1)
        try:
            Chem.SanitizeMol(this_mol)
        except:
            print('Could not sanitize molecule', this_mol.GetProp('_Name'))
        mol_list.append(this_mol)
    return mol_list


#0 RMSD 
def calculate_RMSD(GRO,XTC,ligand_name="MOL"):
    print("Calculating RMSD")
    u = mda.Universe(GRO, XTC, refresh_offsets=True)

    rmsd_array=[]
    
    ligand = u.select_atoms("resname MOL")
    ligand_heavy_atoms= u.select_atoms("resname MOL and not (name H*)")

    protein = u.select_atoms("protein")
    protein_heavy_atoms = u.select_atoms("protein and not (name H*)")
    
    # Calculate 
    RMSD_protein = rms.RMSD(protein_heavy_atoms, protein_heavy_atoms, select='protein',weights='mass',  ref_frame=0)
    RMSD_protein.run()

    RMSD_ligand = rms.RMSD(ligand_heavy_atoms, ligand_heavy_atoms, select='resname MOL',weights='mass', ref_frame=0)
    RMSD_ligand.run()

    df_protein = pd.DataFrame(RMSD_protein.rmsd,columns=['Frame', 'Time (ns)','PTN'])
    df_ligand = pd.DataFrame(RMSD_ligand.rmsd,columns=['Frame', 'Time (ns)','LIG'])

    # Calculate Stats 
    protein_rmsd_avg=(df_protein["PTN"].mean())
    protein_rmsd_std=(df_protein["PTN"].std())
    ligand_rmsd_avg=(df_ligand["LIG"].mean())
    ligand_rmsd_std=(df_ligand["LIG"].std())
    
    rmsd_array.append(protein_rmsd_avg)
    rmsd_array.append(protein_rmsd_std)
    rmsd_array.append(ligand_rmsd_avg)
    rmsd_array.append(ligand_rmsd_std)
    
    #return(protein_rmsd_avg,protein_rmsd_std,ligand_rmsd_avg,ligand_rmsd_std)
    #print(rmsd_array)
    print("Done")
    return rmsd_array


def check_ligand_displacement(GRO,XTC,ligand_name="MOL"):
    
    print("Calculating Ligand center of mass Displacement")

    u = mda.Universe(GRO, XTC, refresh_offsets=True)
    
    ligand = u.select_atoms("resname MOL and not (name H*)")
    ref=(ligand.center_of_mass())

    com_array=[]
    ligand_displacement=[]
    for ts in u.trajectory:
        com_array.append((norm(ref-ligand.center_of_mass())))
    com_array=np.array(com_array)
#
    ligand_displacement.append(com_array.mean())
    ligand_displacement.append(com_array.std())
    
    #print(ligand_displacement)\
    print("Done")
    return ligand_displacement



#1 Process gromacs trajectory

def prepare_inputs():
    file=open("mmbpsa.in","w")
    file.write("&general\n")
    file.write("   verbose=1,\n")
    file.write("   entropy=0,\n")
    file.write("   interval=5,\n")
    file.write("   netcdf=1,\n")
    file.write("   use_sander=1,\n")
    file.write("/\n")
    file.write("&gb\n")
    file.write("  igb=5, saltcon=0.150,\n")
    file.write("/\n")
    file.close()

    file2=open("closest20.in","w")
    file2.write("trajin energy.ncdf\n")
    file2.write("solvent :SOL\n")
    file2.write("closestwaters 20 :MOL oxygen\n")
    file2.write("strip :NA,CL,Na+,Cl- outprefix N20\n")
    file2.write("trajout N20.nc\n")
    file2.write("go\n")
    file2.close()
#Input file for running PB and GB

prepare_inputs()
#sys.exit()

def box_size():
    b = (subprocess.getoutput("tail -n 1 prod.gro").split()[0])
    b = float(b)
    b = int(b)
    return b


def mol_index():

    a=subprocess.getoutput("echo q | /home/magarkar/anaconda3/envs/gromacs/bin/gmx make_ndx -f prod.tpr").split("\n")
    for line in a:
        if "MOL" in line:
            return(line.split()[0])

def postprocess_gromacs(xtc):
    if os.path.isfile(xtc) and os.path.isfile('prod.tpr'):
        mol_num=mol_index()
        b=box_size()*2
        b=str(b)

        subprocess.getoutput("echo %s System| gmx trjconv -f %s -s prod.tpr -center -pbc mol -o center.xtc"%(mol_num,xtc))
#        subprocess.getoutput("echo C-alpha System| gmx trjconv -f center.xtc -s prod.tpr -o traj.xtc -fit rot+trans")
#        subprocess.getoutput("echo C-alpha System| gmx trjconv -f center.xtc -s prod.tpr -o traj.gro -e 0 -fit rot+trans")
        subprocess.getoutput("echo C-alpha System| gmx trjconv -f center.xtc -s prod.tpr -o box.xtc -fit rot+trans")
        subprocess.getoutput("echo Protein System| gmx trjconv -f box.xtc -s prod.tpr -box %s %s %s -o traj.xtc -center"%(b,b,b))
        subprocess.getoutput("echo Protein System| gmx trjconv -f box.xtc -s prod.tpr -box %s %s %s -o traj.gro -e 0 -center"%(b,b,b))
    else:
        print("file not found")

print("Step 1/6 Postprocessing trajectory")

gro='prod.gro'
xtc='prod.xtc'
u=mda.Universe(gro,xtc)

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

#sys.exit()

########################################################################################

#2 parmed to convert gromacs to AMBER
print("Step 2/6 Convering to AMBER topology")

top = pmd.load_file("topol.top",defines=dict(DEFINE='DUMMY'))
struct = pmd.load_file("traj.gro",skip_bonds=True)

top.save('complex.parm7',overwrite=True)
struct.save('complex.rst7',overwrite=True)

#3 MDAnalysis xtc to netcdf
print("Step 3/6 Convering to AMBER trajectory")

u = mda.Universe('traj.gro', 'traj.xtc')
u.atoms.write('energy.ncdf', frames='all')

#4 Prepare trajectory with waters
print("Step 4/6 Preparing trajectory with waters")

subprocess.getoutput("cpptraj -p complex.parm7 -i closest20.in")

#5 Separate components
print("Step 5/6 Separating topology components")

subprocess.getoutput('ante-MMPBSA.py -p N20.complex.parm7 -c complex_nowat.prmtop -r protein.prmtop -l ligand.prmtop -s ":NA,CL,Na+,Cl-" -n ":MOL" --radii=mbondi2') 

#6 Calculate DG
print("Step 6/6 Running MMGBSA")

subprocess.getoutput("MMPBSA.py -O -i mmbpsa.in -o result.dat -sp N20.complex.parm7 -cp complex_nowat.prmtop -rp protein.prmtop -lp ligand.prmtop -y N20.nc")

subprocess.getoutput("rm _MM* *#*")
#subprocess.getoutput("rm whole.xtc nojump.xtc energy.ncdf N20.nc")
u2 = mda.Universe("N20.complex.parm7","N20.nc")
u2.atoms.write('for_simba_raw.xtc', frames='all')
u2.atoms.write('for_simba.gro', frames=[0])

subprocess.getoutput("echo C-alpha System | gmx trjconv -f for_simba_raw.xtc -s for_simba.gro -o for_simba.xtc -fit rot+trans")
subprocess.getoutput("echo C-alpha System | gmx trjconv -f for_simba_raw.xtc -s for_simba.gro -o for_simba.pdb -e 0 -fit rot+trans")

#
ligand_displacement = check_ligand_displacement("for_simba.gro","for_simba.xtc",ligand_name="MOL")
rmsd_summary=calculate_RMSD("for_simba.gro","for_simba.xtc")

subprocess.getoutput("egrep \"ATOM|HETATM\" %s | egrep %s > ligand.pdb"%("for_simba.pdb","MOL"))
subprocess.getoutput("/home/magarkar/anaconda3/bin/obabel -ipdb ligand.pdb -O ligand.sdf --gen2d")
#
mol = load_sdf_with_N("ligand.sdf")
mol = mol[0]
mol = Chem.RemoveHs(mol)
mol.SetProp('Ligand_Displacement',str(ligand_displacement))
mol.SetProp('RMSD_Summary(PTN/stddev/LIG/stddev)',str(rmsd_summary))
if os.path.isfile('result.dat'):
    dg=subprocess.getoutput("egrep 'DELTA TOTAL' result.dat").split()[2]
    mol.SetProp('dg',str(dg))
w = Chem.SDWriter('SIMBA.sdf')
w.write(mol)
w.close()
#

if os.path.isfile("traj.xtc") and os.path.isfile("traj.gro"):
    print("Cleaning UP")
    subprocess.getoutput("rm *#* *.prmtop mmbpsa.in closest20.in complex.* energy.* _MM* N20.* traj.xtc short.xtc index.ndx center.xtc for_simba_raw.xtc")
