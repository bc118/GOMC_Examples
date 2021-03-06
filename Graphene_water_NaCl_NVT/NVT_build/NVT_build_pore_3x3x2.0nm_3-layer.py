from files.porebuilder import GraphenePore
import mbuild as mb
from foyer import Forcefield
import mbuild.formats.charmm_writer as mf_charmm

Water_mol2_file = 'files/tip3p.mol2'
Fake_water_mol2_file = 'files/fake_tip3p.mol2'

Water_res_name = 'H2O'
Fake_water_res_name = 'h2o'

FF_file_water_graphene_NaCl = 'files/FF_graphene_SPCE_NaCl.xml'
FF_file_fake_water = 'files/FF_Fake_SPCE.xml'

water = mb.load(Water_mol2_file)
water.name = Water_res_name
water.energy_minimization(forcefield = FF_file_water_graphene_NaCl  , steps=10**9)

fake_water = mb.load(Fake_water_mol2_file )
fake_water.name = Fake_water_res_name
fake_water.energy_minimization(forcefield = FF_file_fake_water , steps=10**9)


Molecule_Na = mb.Compound(name="Na")
Molecule_Na.name = 'Na'  # naming the molecule (i.e., residue)

Molecule_Cl = mb.Compound(name="Cl")
Molecule_Cl.name = 'Cl'  # naming the molecule (i.e., residue)




FF_Graphene_pore_w_solvent_Dict = {Molecule_Na.name : FF_file_water_graphene_NaCl,
                                   Molecule_Cl.name : FF_file_water_graphene_NaCl,
                                   water.name: FF_file_water_graphene_NaCl ,
                                   fake_water.name: FF_file_fake_water,
                                   'BOT' : FF_file_water_graphene_NaCl ,
                                   'TOP' : FF_file_water_graphene_NaCl
                                  }

residues_Graphene_pore_w_solvent_List = [ Molecule_Na.name,
                                         Molecule_Cl.name,
                                          water.name,
                                         fake_water.name,
                                         'BOT', 'TOP']

Fix_bonds_angles_residues = [ water.name, fake_water.name]

Fix_Graphene_residue = [ 'BOT', 'TOP']



pore_width_nm = 2.0
No_sheets = 3
sheet_spacing = 0.335

water_spacing_from_walls = 0.2

Total_molecules = 500
n_fake_waters = 5
n_Na = 10
n_Cl = 10

#for GOMC, currently we need to add the space at the end of the simulation
# this does not matter as we are using PBC's
empty_graphene_pore = GraphenePore(
        pore_width=sheet_spacing ,
        pore_length=3.0,
        pore_depth=3.0,
        n_sheets=No_sheets,
        slit_pore_dim=2
)

empty_graphene_pore_shifted = empty_graphene_pore


n_waters = Total_molecules - n_fake_waters - n_Na - n_Cl

#note the default spacing of 0.2 automatically accounted for in the water box packing (i.e. adding 0.2 nm for 1 wall is really 0.4 nm)
water_between_pores = mb.fill_box(compound=[water,fake_water, Molecule_Na, Molecule_Cl ],
                                  n_compounds= [n_waters, n_fake_waters, n_Na, n_Cl] ,
                                  box=[3.0, 3.0, pore_width_nm - water_spacing_from_walls*1])
water_between_pores.translate([0,  0, sheet_spacing*(2*No_sheets-1) + water_spacing_from_walls])

filled_pore = empty_graphene_pore
filled_pore.add(water_between_pores, inherit_periodicity=False)
filled_pore.translate([ -filled_pore.center[0],   -filled_pore.center[1], 0])
filled_pore.periodicity[2] = sheet_spacing*(2*No_sheets-1)+pore_width_nm





print('Running: GOMC FF file, and the psf and pdb files')
mf_charmm.charmm_psf_psb_FF(filled_pore,
                            'filled_pore_water_fake_water_NaCl_3x3x2.0nm_3-layer',
                            structure_1 =None ,
                            filename_1 = None,
                            GOMC_FF_filename ="GOMC_pore_water_fake_water_NaCl_FF" ,
                            forcefield_files = FF_Graphene_pore_w_solvent_Dict,
                            residues= residues_Graphene_pore_w_solvent_List ,
                            Bead_to_atom_name_dict = None,
                            fix_residue = Fix_Graphene_residue,
                            fix_res_bonds_angles = Fix_bonds_angles_residues,
                            reorder_res_in_pdb_psf = False
                            )
print('Completed: GOMC FF file, and the psf and pdb files')