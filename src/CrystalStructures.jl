module CrystalStructures

using   Printf
using   StaticArrays
using   LinearAlgebra
import  Crystalline

# Number of space groups
# OEIS: A006227
const NO_SGS = (2, 17, 230, 4985, 222097)
# Number of space groups, with chiral copies not distinct
# OEIS: A004029
const NO_SGS_PAIRED = (2, 17, 219, 4783, 222018, 28927915)

# List of elements
# You can convert a vector of element numbers to a vector of element names with ELEMENTS[v]
const ELEMENTS = 
[ 
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "As", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
]

# TODO: can we modify this to make it unnecessary?
const ELEMENT_LOOKUP = Dict{String, Int}([(ELEMENTS[n] => n) for n in 1:length(ELEMENTS)])

# Metallic radii of the elements in angstroms
# Probably will be useful for sphere packing related things
const RADII_METALLIC_ANGSTROMS = Dict{String, Float64}(
    "Li" => 1.52,
    "Be" => 1.12,
    "Na" => 1.86,
    "Mg" => 1.60,
    "Al" => 1.43,
    "K"  => 2.27,
    "Ca" => 1.97,
    "Sc" => 1.62,
    "Ti" => 1.47,
    "V"  => 1.34,
    "Cr" => 1.28,
    "Mn" => 1.27,
    "Fe" => 1.26,
    "Co" => 1.25,
    "Ni" => 1.24,
    "Cu" => 1.28,
    "Zn" => 1.34,
    "Ga" => 1.35,
    "Rb" => 2.48,
    "Sr" => 2.15,
    "Y"  => 1.80,
    "Zr" => 1.60,
    "Nb" => 1.46,
    "Mo" => 1.39,
    "Tc" => 1.36,
)

include("vectors.jl")
include("datagrids.jl")
include("crystals.jl")
include("filetypes.jl")

# Types
export AbstractCrystal, CrystalStructure, AbstractCrystalData, CrystalStructureWithData
# Functions
export formula, readXYZ, writeXYZ, readXSF3D, writeXSF3D

end # module