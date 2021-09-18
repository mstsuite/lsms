systemid = "fe"
system_title = "Iron with LSF"
pot_in_type = 1
pot_out_type = 1

num_atoms = 1
nspin = 2
mtasa = 0
temperature = 300

xcFunctional = { 0, 1 } -- von Barth Hedin (LSMS_1)

iprint = -1
default_iprint = -1
print_node = 0
istop = "main"

nscf = 50
rmsTolerance = 1.0e-9

xRepeat = 1
yRepeat = 1
zRepeat = 1
makeTypesUnique = 1

energyContour = { npts = 25, grid = 2, ebot = -0.5, etop = 0.0, eitop = 0.825, eibot = 0.0025 }

a = 6.75

bravais = {}
bravais[1] = { 0.5 * a, 0.5 * a, 0 }
bravais[2] = { 0, 0.5 * a, 0.5 * a }
bravais[3] = { 0.5 * a, 0, 0.5 * a }

site_default = { lmax = 3,
                 rLIZ = 8.5,
                 rsteps = { 89.5, 91.5, 93.2, 99.9 },
                 atom = "Fe",
                 Z = 26, Zc = 10, Zs = 8, Zv = 8, rad = 2, lsf = 0 }

mixing = { { quantity = "potential", algorithm = "broyden", mixing_parameter = 0.05 } }

numberOfMixQuantities = 0
for k, v in pairs(mixing) do
    numberOfMixQuantities = numberOfMixQuantities + 1
end

site = {}
for i = 1, num_atoms do
    site[i] = {}
end

site[1].pos = { 0, 0, 0 }
site[1].evec = { 0, 0, 1 }

-- set site defaults
for i = 1, num_atoms do
    for k, v in pairs(site_default) do
        if (site[i][k] == nil) then
            site[i][k] = v
        end
    end
end
