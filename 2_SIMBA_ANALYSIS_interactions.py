#!/home/magarkar/anaconda3/envs/simba/bin/python

import warnings
warnings.filterwarnings("ignore")
# Basics
import subprocess,os,operator,sys
from tqdm import tqdm_notebook
from tqdm import tqdm

#Statistics
from collections import Counter
import numpy as np

from plip.structure.preparation import PDBComplex
from plip.exchange.report import BindingSiteReport
from plip.basic import config
config.NOHYDRO = True
config.SILENT = True
config.QUIET = True
config.VERBOSE = False

import pandas as pd

#MD
import MDAnalysis as mda
from MDAnalysis.analysis import rms

#RDKIT
from rdkit import Chem,Geometry
from rdkit.Chem import Draw,rdDepictor
from rdkit.Chem.Draw import SimilarityMaps
import matplotlib
matplotlib.use('Agg')


#PDB in DF
from biopandas.pdb import PandasPdb

def generate_rmsf_image():
    mol=Chem.MolFromMolFile('mol.sdf')
    ligand_rmsf_arr=calculate_ligand_RMSF(pdb_file,xtc_file)
    min_value = np.amin (ligand_rmsf_arr)
    ligand_rmsf_arr = ligand_rmsf_arr - min_value
    max_value = np.amax (ligand_rmsf_arr)
    ligand_rmsf_arr = ligand_rmsf_arr / max_value
    ligand_rmsf_arr = ligand_rmsf_arr - 0.5
    
    rdDepictor.Compute2DCoords(mol)

    img=SimilarityMaps.GetSimilarityMapFromWeights(mol,ligand_rmsf_arr.tolist(),alpha=0.2,colorMap="coolwarm",size=(600,600));
    img.savefig("ligand_rmsf.png",bbox_inches="tight")    

def calculate_ligand_RMSF(GRO,XTC,ligand_name="MOL"):
    print("Calculating Ligand RMSF")
    u = mda.Universe(GRO, XTC,refresh_offsets=True)
    selected_atoms = u.select_atoms("resname MOL and not (name H*)")
    rmsfer = rms.RMSF(selected_atoms).run()
    ligand_rmsf_array=rmsfer.rmsf
    print("Done")
    return ligand_rmsf_array

def simba_figures(polar_summary_num,polar_summary_res,nonpolar_summary_num,nonpolar_summary_res):
    mol=Chem.MolFromMolFile('mol.sdf')
    polar_mol=mol
    rdDepictor.Compute2DCoords(polar_mol)

    for i in range(len(polar_summary_res)):
        if polar_summary_res[i] != "":
            polar_mol.GetAtomWithIdx(i).SetProp("atomNote",polar_summary_res[i])

    obj = my_draw(polar_mol,polar_summary_num, size= (600,600))
    obj.FinishDrawing()
    obj.WriteDrawingText("polar.png")

    mol=Chem.MolFromMolFile('mol.sdf')
    nonpolar_mol=mol

    rdDepictor.Compute2DCoords(nonpolar_mol)

    for i in range(len(nonpolar_summary_res)):
        if nonpolar_summary_res[i] != "":
            nonpolar_mol.GetAtomWithIdx(i).SetProp("atomNote",nonpolar_summary_res[i])

    obj = my_draw(nonpolar_mol,nonpolar_summary_num, size= (600,600))
    obj.FinishDrawing()
    obj.WriteDrawingText("nonpolar.png")

def my_draw(mol, weights, colorMap=None, scale=-1, size=(600, 600),
                                sigma=None, coordScale=1.5, step=0.001, colors='k', contourLines=5,
                                alpha=0.2, draw2d=None, **kwargs):
    if mol.GetNumAtoms() < 2:
        raise ValueError("too few atoms")
    
    draw2d = Draw.MolDraw2DCairo(size[0], size[1])
    draw2d.SetFontSize(0.5)
    draw2d.SetLineWidth(1)
#    draw2d.SetLineWidth(2)
    
    if sigma is None:
        if mol.GetNumBonds() > 0:
            bond = mol.GetBondWithIdx(0)
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            sigma = 0.3 * (mol.GetConformer().GetAtomPosition(idx1)-mol.GetConformer().GetAtomPosition(idx2)).Length()
        else:
            sigma = 0.3 * (mol.GetConformer().GetAtomPosition(0)-mol.GetConformer().GetAtomPosition(1)).Length()
        sigma = round(sigma, 2)
    sigmas = [sigma]*mol.GetNumAtoms()
    locs=[]
    for i in range(mol.GetNumAtoms()):
        p = mol.GetConformer().GetAtomPosition(i)
        locs.append(Geometry.Point2D(p.x,p.y))
    draw2d.ClearDrawing()
    ps = Draw.ContourParams()
    ps.fillGrid=True
    ps.gridResolution=0.1
    ps.extraGridPadding = 0.5
    Draw.ContourAndDrawGaussians(draw2d,locs,weights,sigmas,nContours=contourLines,params=ps)

    draw2d.drawOptions().clearBackground = False
    draw2d.DrawMolecule(mol)
    return draw2d

def generate_sdf(pdb_file):
    input_pdb=open(pdb_file,"r")
    mol_pdb=open('mol.pdb',"w")
    for line in input_pdb.readlines():
        if 'MOL' in line and 'ATOM' in line:
            mol_pdb.write(line)
    mol_pdb.close()
    subprocess.getoutput("/home/magarkar/anaconda3/bin/obabel -ipdb mol.pdb -O mol.sdf -d")

# Post process the DF and generate a summary

def post_porcess_df(df):
    df_num=(~df.isna()).astype(int)
    df_num_list=(df_num.mean().to_list())
    df_num_list = [ round(elem, 1) for elem in df_num_list ]

    df_str=df.fillna("")
    intertation_str=[]
    
    for column in df_str.columns:
        column_str=""
        temp_list=(df_str[column].to_list())
        temp_str=(','.join(temp_list))
        temp_list2 = temp_str.split(",")
        temp_list2 = [i for i in temp_list2 if i]
        frequency = Counter(temp_list2)
        temp_dict=(dict(frequency))
        sorted_d = dict( sorted(temp_dict.items(), key=operator.itemgetter(1),reverse=True))
        
        max_count=0
        
        for key, value in sorted_d.items():
            max_count+=1
            if max_count<3:
                column_str+=key
                column_str+=","
        
        column_str = column_str[:-1]
        intertation_str.append(column_str)
    
    return df_num_list,intertation_str

# Process the inputs
    
def process_trajectory(gro_file,xtc_file,polar_interactions,nonpolar_interactions):
    u=mda.Universe(gro_file,xtc_file)
    trajectory_len=u.trajectory[-1].time
    print("trajectory len =",trajectory_len/1000,"ns")
    a=u.select_atoms('all')
    #for i in range(0,u.trajectory.n_frames):
    for i in tqdm(range(0,u.trajectory.n_frames)):
        if(os.path.isfile('temp_frame.pdb')):
            subprocess.getoutput('rm temp_frame.pdb')

        a.write('temp_frame.pdb',frames=[i])
        polar_df,nonpolar_df = prepare_the_plip_fp("temp_frame.pdb")
        polar_interactions=polar_interactions.append(polar_df,ignore_index=True)
        nonpolar_interactions=nonpolar_interactions.append(nonpolar_df,ignore_index=True)
        #print(i)
    return polar_interactions,nonpolar_interactions


def prepare_the_plip_fp(input_pdb_file):
    interactions_by_site = retrieve_plip_interactions(input_pdb_file)
    selected_site_list = list(interactions_by_site.keys())
    selected_site=''
    
    #print(selected_site_list)

    for i in range(0,len(selected_site_list)):
        if 'MOL' in selected_site_list[i]:
            #print(i)
            #print(selected_site_list[i])
            selected_site=i
            break


    selected_site=selected_site_list[selected_site]
    #print(selected_site)

    hydrophobic=create_df_from_binding_site(interactions_by_site[selected_site], interaction_type="hydrophobic")
    pistacking=create_df_from_binding_site(interactions_by_site[selected_site], interaction_type="pistacking")

    hydrophobic.to_html("before.html")
    pistacking.to_html("before_stack.html")

    hbond=create_df_from_binding_site(interactions_by_site[selected_site], interaction_type="hbond")
    waterbridge=create_df_from_binding_site(interactions_by_site[selected_site], interaction_type="waterbridge")
    saltbridge=create_df_from_binding_site(interactions_by_site[selected_site], interaction_type="saltbridge")
    pication=create_df_from_binding_site(interactions_by_site[selected_site], interaction_type="pication")
    halogen=create_df_from_binding_site(interactions_by_site[selected_site], interaction_type="halogen")

    metal=create_df_from_binding_site(interactions_by_site[selected_site], interaction_type="metal")

    #hbond['LIGINFO']=np.where(hbond['PROTISDON']=='True',hbond.ACCEPTORIDX,hbond.DONORIDX)
    #print(hbond)
    hbond['LIGINFO'] = hbond['ACCEPTORIDX'].where(hbond['PROTISDON'].astype(str)=='True',other=hbond['DONORIDX'])
    waterbridge['LIGINFO'] = waterbridge['ACCEPTOR_IDX'].where(waterbridge['PROTISDON'].astype(str)=='True',other=waterbridge['DONOR_IDX'])
    halogen['LIGINFO']=halogen['DON_IDX']
    #----
    hbond=post_porcess_simple_interaction_df(hbond)
    waterbridge=post_porcess_simple_interaction_df(waterbridge)
    halogen=post_porcess_simple_interaction_df(halogen)
    saltbridge = post_porcess_complex_interaction_df(saltbridge)
    pication = post_porcess_complex_interaction_df(pication)

    hydrophobic['LIGINFO']=hydrophobic.LIGCARBONIDX
    hydrophobic=post_porcess_simple_interaction_df(hydrophobic)
    pistacking=post_porcess_complex_interaction_df(pistacking)

    hydrophobic.to_html("phobe.html")
    pistacking.to_html("stack.html")
    


    #---------------------------------------------------------
    polar_df=pd.DataFrame(columns=['LIGINFO','INFO'])
    polar_df=polar_df.append(hbond,ignore_index=True)
    polar_df=polar_df.append(halogen,ignore_index=True)
    polar_df=polar_df.append(saltbridge,ignore_index=True)
    polar_df=polar_df.append(waterbridge,ignore_index=True)
    polar_df=polar_df.append(pication,ignore_index=True)

    polar_df.LIGINFO=polar_df.LIGINFO.astype(int)

    nonpolar_df=pd.DataFrame(columns=['LIGINFO','INFO'])
    
    
    nonpolar_df=nonpolar_df.append(hydrophobic,ignore_index=True)
    nonpolar_df=nonpolar_df.append(pistacking,ignore_index=True)

    polar_df.LIGINFO=polar_df.LIGINFO.astype(int)
    nonpolar_df.LIGINFO=nonpolar_df.LIGINFO.astype(int)

    polar_df=polar_df.drop_duplicates()
    nonpolar_df=nonpolar_df.drop_duplicates()

    polar_df=polar_df.reset_index()
    nonpolar_df=nonpolar_df.reset_index()

    polar_df=polar_df.groupby('LIGINFO')['INFO'].apply(','.join).reset_index()
    nonpolar_df=nonpolar_df.groupby('LIGINFO')['INFO'].apply(','.join).reset_index()

    polar_df_dict=polar_df.set_index(polar_df.LIGINFO).to_dict()['INFO']
    nonpolar_df_dict=nonpolar_df.set_index(nonpolar_df.LIGINFO).to_dict()['INFO']
    #---------------------------------------------------------
    return polar_df_dict,nonpolar_df_dict

def create_df_from_binding_site(selected_site_interactions, interaction_type="hbond"):
    """
    Creates a data frame from a binding site and interaction type.

    Parameters
    ----------
    selected_site_interactions : dict
        Precaluclated interactions from PLIP for the selected site
    interaction_type : str
        The interaction type of interest (default set to hydrogen bond).

    Returns
    -------
    DataFrame :
        Data frame with information retrieved from PLIP.
    """

    # check if interaction type is valid:
    valid_types = [
        "hydrophobic",
        "hbond",
        "waterbridge",
        "saltbridge",
        "pistacking",
        "pication",
        "halogen",
        "metal",
    ]

    if interaction_type not in valid_types:
        print("!!! Wrong interaction type specified. Hbond is chosen by default!!!\n")
        interaction_type = "hbond"

    df = pd.DataFrame.from_records(
        # data is stored AFTER the column names
        selected_site_interactions[interaction_type][1:],
        # column names are always the first element
        columns=selected_site_interactions[interaction_type][0],
    )
    
    #
    df["INFO"]=df["RESNR"].astype(str)+df["RESTYPE"].astype(str)

    #
    
    return df

def post_porcess_simple_interaction_df(df_in):
    df_out=pd.DataFrame(columns=['LIGINFO','INFO'])

    for i in range(0,len(df_in)):
        df_out=df_out.append({'LIGINFO':df_in.LIGINFO[i],'INFO':df_in.INFO[i]},ignore_index=True)
    
    df_out=df_out.drop_duplicates()
    return df_out

def post_porcess_complex_interaction_df(df_in):
    df_out=pd.DataFrame(columns=['LIGINFO','INFO'])
    for i in range(0,len(df_in)):
        if "," in df_in.LIG_IDX_LIST[i]:
            lig_atoms=df_in.LIG_IDX_LIST[i].split(",")
            for lig_atom in lig_atoms:
                #print(lig_atom,pistacking.INFO[i])
                df_out=df_out.append({'LIGINFO': lig_atom,'INFO':df_in.INFO[i]}, ignore_index=True)
        else:
            df_out=df_out.append({'LIGINFO': df_in.LIG_IDX_LIST[0],'INFO':df_in.INFO[i]}, ignore_index=True)
            
    df_out=df_out.drop_duplicates()
    return df_out
def retrieve_plip_interactions(pdb_file):
    """
    Retrieves the interactions from PLIP.

    Parameters
    ----------
    pdb_file :
            The PDB file of the complex.

    Returns
    -------
    dict :
            A dictionary of the binding sites and the interactions.
    """
    protlig = PDBComplex()
    protlig.load_pdb(pdb_file)  # load the pdb file
    for ligand in protlig.ligands:
        protlig.characterize_complex(ligand)  # find ligands and analyze interactions
    sites = {}
    # loop over binding sites
    for key, site in sorted(protlig.interaction_sets.items()):
        binding_site = BindingSiteReport(site)  # collect data about interactions
        # tuples of *_features and *_info will be converted to pandas data frame
        keys = (
            "hydrophobic",
            "hbond",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
        )
        # interactions is a dictionary which contains relevant information for each
        # of the possible interactions: hydrophobic, hbond, etc. in the considered
        # binding site. Each interaction contains a list with
        # 1. the features of that interaction, e.g. for hydrophobic:
        # ('RESNR', 'RESTYPE', ..., 'LIGCOO', 'PROTCOO')
        # 2. information for each of these features, e.g. for hydrophobic
        # (residue nb, residue type,..., ligand atom 3D coord., protein atom 3D coord.)
        interactions = {
            k: [getattr(binding_site, k + "_features")] + getattr(binding_site, k + "_info")
            for k in keys
        }
        sites[key] = interactions
    return sites


def prepare_ligand_dataframe(pdb_file):
    pandas_pdb=PandasPdb()
    pandas_pdb=pandas_pdb.read_pdb(pdb_file)


    ligand_df_all = pandas_pdb.df['ATOM'][pandas_pdb.df['ATOM']['residue_name'] == 'MOL']
    #ligand_df_all['atom_number','atom_name']
    ligand_df=ligand_df_all[['atom_number','atom_name']]
    ligand_df=ligand_df[~ligand_df.atom_name.str.contains('H')]
    return ligand_df

###############################################################################################

pdb_file="for_simba.pdb"
xtc_file="for_simba.xtc"

ligand_df=prepare_ligand_dataframe(pdb_file)
MOL_atom_indices=ligand_df['atom_number'].to_list()
generate_sdf(pdb_file)

polar_interactions = pd.DataFrame(columns=MOL_atom_indices)
nonpolar_interactions = pd.DataFrame(columns=MOL_atom_indices)

polar_interactions_df,nonpolar_interactions_df=process_trajectory(pdb_file,xtc_file,polar_interactions,nonpolar_interactions)

polar_summary_num,polar_summary_res=post_porcess_df(polar_interactions_df)
nonpolar_summary_num,nonpolar_summary_res=post_porcess_df(nonpolar_interactions_df)

simba_figures(polar_summary_num,polar_summary_res,nonpolar_summary_num,nonpolar_summary_res)
generate_rmsf_image()
