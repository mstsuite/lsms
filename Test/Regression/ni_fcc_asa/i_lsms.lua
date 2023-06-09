--[[

Test for non-spin polarized TiAl

]]--

systemid = "al"
system_title = "TiAl test for LSMS 3"
pot_in_type = 1
pot_out_type = 1

relativity = "s"
core_relativity = "f"

num_atoms = 2
nspin = 1
mtasa = 1
use_voronoi = 0

xcFunctional = { 0, 2, 2 } -- Vosko-Wilk-Nusair

iprint = -1
default_iprint = -2
print_node = 0
istop = "main"

nscf = 10
rmsTolerance = 1.0e-16

xRepeat = 1
yRepeat = 1
zRepeat = 1
makeTypesUnique = 1

energyContour = { npts = 34, grid = 2, ebot = -0.3, etop = 0.0, eitop = 0.825, eibot = 0.0025 }

lat = 3.48
a = lat / 0.529177

bravais = {}
bravais[1] = { 1.0 * a, 0.0 * a, 0.0 * a }
bravais[2] = { 0.0 * a, 1.0 * a, 0.0 * a }
bravais[3] = { 0.0 * a, 0.0 * a, 1.0 * a }

site_default = { lmax = 3, rLIZ = 6.5,
                 rsteps = { 89.5, 91.5, 93.2, 99.9 },
                 rad = 2, evec = { 0, 0, 1} }

mixing = { 
        { quantity = "potential", algorithm = "simple", alpha = 0.5 },
        { quantity = "density", algorithm = "simple", alpha = 0.01 },
    }


site = {}
for i = 1, num_atoms do
    site[i] = {}
end

site[1].atom = "Ni"
site[1].Z = 28
site[1].Zc = 10
site[1].Zs = 10
site[1].Zv = 8
site[1].pos = { 0.5 * a, 0.5 * a, 0.5 * a }

site[2].atom = "Ni"
site[2].Z = 28
site[2].Zc = 10
site[2].Zs = 10
site[2].Zv = 8
site[2].pos = { 0.5 * a, 0.5 * a, 0.5 * a }

for i = 1, num_atoms do
    for k, v in pairs(site_default) do
        if (site[i][k] == nil) then
            site[i][k] = v
        end
    end
end

