--[[

Test for non-spinpolarized Iron in a conventional bcc cell

]]--

systemid = "fe"
system_title = "Iron test for LSMS 3"
pot_in_type = 1
pot_out_type = 1

relativity = "s"
core_relativity = "f"

num_atoms = 2
nspin = 1
mtasa = 0

xcFunctional = { 1, 1, 7 }

iprint = -1
default_iprint = -1
print_node = 0
istop = "main"

nscf = 5
rmsTolerance = 1.0e-16

xRepeat = 1
yRepeat = 1
zRepeat = 1
makeTypesUnique = 1

energyContour = { npts = 31, grid = 2, ebot = -0.3, etop = 0.0, eitop = 0.825, eibot = 0.0025 }

a = 3.00 / 0.529177

bravais = {}
bravais[1] = { 1.00 * a, 0.01 * a, 0.01 * a }
bravais[2] = { 0.00 * a, 0.98 * a, 0.00 * a }
bravais[3] = { 0.00 * a, 0.00 * a, 1.05 * a }

site_default = { lmax = 2, rLIZ = 4.5, rsteps = { 89.5, 91.5, 93.2, 99.9 },
                 atom = "Fe", Z = 26, Zc = 10, Zs = 8, Zv = 8, rad = 2 }

mixing = { { quantity = "potential", algorithm = "simple", mixing_parameter = 0.05 } }

numberOfMixQuantities = 0
for k, v in pairs(mixing) do
    numberOfMixQuantities = numberOfMixQuantities + 1
end

site = {}
for i = 1, num_atoms do
    site[i] = {}
end


-- site 1: Fe
site[1].pos = { 0, 0, 0.05 * a }
site[1].evec={0,0,1}
site[1].pot_in_idx=0
site[1].atom="Fe"
site[1].Z=26
site[1].Zc=10
site[1].Zs=8
site[1].Zv=8
-- site 2: Pt
site[2].pos = { 0.45 * a, 0.52 * a, 0.48 * a }
site[2].evec={0,0,1}
site[2].pot_in_idx=1
site[2].atom="Pt"
site[2].Z=78
site[2].Zc=46
site[2].Zs=22
site[2].Zv=10

for i = 1, num_atoms do
    for k, v in pairs(site_default) do
        if (site[i][k] == nil) then
            site[i][k] = v
        end
    end
end


