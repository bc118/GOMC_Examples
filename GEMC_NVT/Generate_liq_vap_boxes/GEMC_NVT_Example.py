import mbuild as mb
from foyer import Forcefield
import mbuild.formats.charmm_writer as mf_charmm


forcefield_names = 'trappe-ua'
FF =  Forcefield(name=forcefield_names)

Molecule_A =mb.load('files/iso_decane.mol2')
FF.apply(Molecule_A)
Molecule_A.name = 'IDE'  # naming the molecule (i.e., residue)

Molecule_B =mb.load('files/iso_octane.mol2')
FF.apply(Molecule_B)
Molecule_B.name = 'IOT'  # naming the molecule (i.e., residue)

approx_total_No_atoms_liq = 10000
approx_total_No_atoms_per_cg_atom = 3.25    # 1 for AA and ~3.0 to 3.5 for UA.
est_beads_per_moledule = 9
est_total_molecules_liq = approx_total_No_atoms_liq / approx_total_No_atoms_per_cg_atom/est_beads_per_moledule

min_atom_spacing = 0.25
vol_vap_div_vol_liq = 10
molecules_vap_div_molecules_liq = 0.1

mol_fraction_Molecule_A = 0.25
mol_fraction_Molecule_B = 0.75

Molecule_Type_List = [Molecule_A, Molecule_B ]
Molecule_ResName_List = [Molecule_A.name, Molecule_B.name ]
Total_molecules_liquid = [int(mol_fraction_Molecule_A * est_total_molecules_liq) ,
                          int(mol_fraction_Molecule_B * est_total_molecules_liq)]
Total_molecules_vapor = [int(mol_fraction_Molecule_A * est_total_molecules_liq/10) ,
                         int(mol_fraction_Molecule_B * est_total_molecules_liq/10)]
Bead_to_atom_name_dict = { '_CH3':'C', '_CH2':'C',  '_CH':'C', '_HC':'C'}

print('Running: liquid phase box packing')
box_liq = mb.fill_box(compound=Molecule_Type_List,
                      n_compounds=Total_molecules_liquid,
                      box=[4, 4, 4])
print('Completed: liquid phase box packing')

print('Running: vapor phase box packing')
box_vap = mb.fill_box(compound=Molecule_Type_List,
                      n_compounds=Total_molecules_vapor,
                      box=[8, 8, 8])
print('Completed: vapor phase box packing')

print('Running: GOMC FF file, and the psf and pdb files')
mf_charmm.charmm_psf_psb_FF(box_liq,
                            'IOT_IDE_liquid_box',
                            structure_1 =box_vap ,
                            filename_1 = 'IOT_IDE_vapor_box',
                            GOMC_FF_filename ="GOMC_IOT_IDE_FF.inp" ,
                            forcefield_names= forcefield_names ,
                            residues= Molecule_ResName_List ,
                            Bead_to_atom_name_dict = Bead_to_atom_name_dict,
                            fix_residue = None,
                            fix_res_bonds_angles = None,
                            reorder_res_in_pdb_psf = False)
print('Completed: GOMC FF file, and the psf and pdb files')
print('Completed: liquid phase box packing')


