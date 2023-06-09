--[[

Test for non-spin polarized Al

]] --
systemid = "al"
system_title = "Al test for LSMS 3"
pot_in_type = -1
pot_out_type = 1

relativity = "s"
core_relativity = "f"

num_atoms = 4
nspin = 1
use_voronoi = 0

xcFunctional = {1, 1, 7} 

iprint = -1
default_iprint = -2
print_node = 0
istop = "main"

efermi = 0.65
rmax = 35.0
h_step = 0.025

nscf = 100
rmsTolerance = 1.0e-9
energyTolerance = 1.0e-8

xRepeat = 1
yRepeat = 1
zRepeat = 1
makeTypesUnique = 1


energyContour = { npts = 24, grid = 2, ebot = -0.3, etop = 0.0, eitop = 0.825, eibot = 0.0025 }

lat = 4.03
a = lat / 0.529177

bravais = {}
bravais[1] = { 1.0 * a, 0.0 * a, 0.0 * a }
bravais[2] = { 0.0 * a, 1.0 * a, 0.0 * a }
bravais[3] = { 0.0 * a, 0.0 * a, 1.0 * a }

site_default = {
    lmax = 3,
    rLIZ = 9.0 * lat / 4.05,
    rsteps = {89.5, 91.5, 93.2, 99.9},
    atom = "Al",
    rad = 2,
    evec = {0, 0, 1}
}

n_init_iterations = 0

debug_chem_pot = false
debug_charge = false
debug_radial_charge = false
debug_madelung = false
debug_potential = false

mixing = {
    {quantity = "potential", algorithm = "broyden", alpha = 0.05}
}

site = {}
for i = 1, num_atoms do
    site[i] = {}
end

site[1].atom = "Al"
site[1].Z = 13
site[1].Zc = 2
site[1].Zs = 8
site[1].Zv = 3

site[2].atom = "Al"
site[2].Z = 13
site[2].Zc = 2
site[2].Zs = 8
site[2].Zv = 3


site[3].atom = "Al"
site[3].Z = 13
site[3].Zc = 2
site[3].Zs = 8
site[3].Zv = 3

site[4].atom = "Al"
site[4].Z = 13
site[4].Zc = 2
site[4].Zs = 8
site[4].Zv = 3

site[1].pos = { -0.25 * a, -0.25 * a, -0.25 * a }
site[2].pos = { 0.25 * a, 0.25 * a, -0.25 * a }
site[3].pos = { -0.25 * a, 0.25 * a, 0.25 * a }
site[4].pos = { 0.25 * a, -0.25 * a, 0.25 * a }

for i = 1, num_atoms do
    for k, v in pairs(site_default) do
        if (site[i][k] == nil) then
            site[i][k] = v
        end
    end
end

