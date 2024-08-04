import unyt as u
from mosdef_dihedral_fit.dihedral_fit.fit_dihedral_with_gomc import fit_dihedral_with_gomc
import os

# The MoSDeF force field (FF) XML file which will be used.
forcefield_file = '../src/ff_xml/oplsaa_ethane_HC_CT_CT_HC.xml'

# The mol2 file which is used.
mol2_file = '../src/starting_coords/ethane_aa.mol2'

# The temperature of the Molecular Mechanics (MM) simulation.
temperature_in_unyt_units = 298.15 * u.Kelvin

# Override the chosen mixing rule (combining_rule), 'geometric' or 'lorentz'),
# which was set in the force field (FF) XML file.  This variable is not required and will be
# selected automatically; however, you can override it if you are unsure of the setting.
combining_rule = 'geometric'

# Atom type naming convention ( str, optional, default=’all_unique’, (‘general’ or ‘all_unique’) )
# Using 'general' is preferred as it makes it easier to create the force field (FF) and makes it 
# and the other generated files more readable.  However, to the 'general' option, the user needs 
# to make sure that all the FF atomclasses in the FF XML file (i.e., atomclass="CT") have the same 
# sigma and epsilon values. Otherwise, it will cause errors.
atom_type_naming_style = 'general'

# The GOMC binary path.
gomc_binary_directory = "path_to_GOMC_folder_GOMC_2_75/bin"

# Load 1 or more Gaussian files (keys), and its value ([0]), which is a list of Gaussian points to remove
# from the fitting process, where the first minimized Gaussian point is removed (i.e., ([0])).
# More than 1 Gaussian file can be loaded, allowing the user to run multiple dihedral angles in separate file,
# minimizing the time required to run the simulations
# (i.e., the user can split them up into many simulations to obtain the full dihedral rotation).
qm_log_file_dict = {
    '../src/guassian_runs_multiplicity_1/HC_CT_CT_HC_multiplicity_1.log': [0],
}

# The dihedral which is being fit.
fit_dihedral_atom_types = ['HC', 'CT', 'CT', 'HC']

# All the other dihedrals which can be zeroed in the fitting process, in a nested list
zero_dihedral_atom_types = None

# Run the "fit_dihedral_with_gomc" command.
fit_dihedral_with_gomc(
    fit_dihedral_atom_types,
    mol2_file,
    forcefield_file,
    temperature_in_unyt_units,
    gomc_binary_directory,
    qm_log_file_dict,
    zero_dihedral_atom_types=zero_dihedral_atom_types,
    qm_engine="gaussian",
    combining_rule=combining_rule,
    atom_type_naming_style=atom_type_naming_style,
    gomc_cpu_cores=1,
    r_squared_min=0.99,
    r_squared_rtol=1e-03
)

# The OPLS dihedral fit constants
print('The OPLS dihedral fit constants:\n')
os.system("cat opls_dihedral_k_constants_fit_energy.txt")

# The converted OPLS to periodic dihedral fit constants
print('The converted OPLS to periodic dihedral fit constants:\n')
os.system("cat periodic_dihedral_k_constants_fit_energy.txt")

# The converted OPLS to RB-torsions dihedral fit constants
print('The converted OPLS to RB-torsions dihedral fit constants:\n')
os.system("cat RB_torsion_k_constants_fit_energy.txt")

# This file contains the raw points used in the fits (QM - MM_less_dihedral(s)_energy being fit)
print('This file contains the raw points used in the fits (QM - MM_less_dihedral(s)_energy being fit):\n')
os.system("cat all_normalized_energies_in_kcal_per_mol.txt")