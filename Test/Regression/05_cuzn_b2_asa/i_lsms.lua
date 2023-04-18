--[[

Test for spin polarized Fe

]] --
systemid = "cuzn"
system_title = "CuZn test for LSMS 3"
pot_in_type = -1
pot_out_type = 1

relativity = "s"
core_relativity = "f"

num_atoms = 2
nspin = 1
mtasa = 1
use_voronoi = 0

xcFunctional = {1, 1, 7}

iprint = -1
default_iprint = -1
print_node = 0
istop = "main"

nscf = 50
rmsTolerance = 1.0e-11
energyTolerance = 1.0e-8

efermi = 0.264552429606
rmax = 35.0
h_step = 0.015
rmin = 0.0001

xRepeat = 1
yRepeat = 1
zRepeat = 1
makeTypesUnique = 1

energyContour = { npts = 28, grid = 2, ebot = -0.3, etop = 0.0, eitop = 0.825, eibot = 0.0025 }

lat = 3.61

a = lat / 0.529177

bravais = {}
bravais[1] = { 1.0 * a, 0.0 * a, 0.0 * a }
bravais[2] = { 0.0 * a, 1.0 * a, 0.0 * a }
bravais[3] = { 0.0 * a, 0.0 * a, 1.0 * a }

site_default = {
    lmax = 3,
    rLIZ = 11.0,
    rsteps = {89.5, 91.5, 93.2, 99.9},
    atom = "Fe",
    rad = 2,
    mag_mom = 2.0,
    evec = {0, 0, 1}
}

debug_atomic = false
debug_madelung = false
debug_core_states = false
debug_radial_charge = false
debug_charge = false
debug_potential = false
debug_energy = false
debug_convergence = false

mixing = {
    {quantity = "potential", algorithm = "broyden", alpha = 0.05},
}

numberOfMixQuantities = 0
for k, v in pairs(mixing) do
    numberOfMixQuantities = numberOfMixQuantities + 1
end

site = {}
for i = 1, num_atoms do
    site[i] = {}
end

site[1].atom = "Cu"
site[1].Z = 29
site[1].Zc = 10
site[1].Zs = 8
site[1].Zv = 11
site[1].pos = { 0.0 * a, 0.0 * a, 0.0 * a }

site[2].atom = "Zn"
site[2].Z = 30
site[2].Zc = 10
site[2].Zs = 8
site[2].Zv = 12
site[2].pos = { 0.5 * a, 0.5 * a, 0.5 * a }

for i = 1, num_atoms do
    for k, v in pairs(site_default) do
        if (site[i][k] == nil) then
            site[i][k] = v
        end
    end
end

