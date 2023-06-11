# Simple rigid body system

log log_${T} append
units           lj
atom_style      full

pair_style      lj/cut 1.4
bond_style harmonic
angle_style harmonic
dihedral_style harmonic

read_data       last_conf


pair_coeff * * 0        0.359  1.4

pair_coeff 1 1 70       0.900  1.00
pair_coeff 2 5 1600     0.359  0.56
pair_coeff 3 6 1600     0.359  0.56
pair_coeff 4 7 1600     0.359  0.56

pair_coeff 8  11 14    0.359  0.56
pair_coeff 9  12 14    0.359  0.56
pair_coeff 10 13 14    0.359  0.56

pair_coeff 8  8  7    0.358  0.56
pair_coeff 9  9  7    0.358  0.56
pair_coeff 10 10 7    0.358  0.56
pair_coeff 11 11 7    0.358  0.56
pair_coeff 12 12 7    0.358  0.56
pair_coeff 13 13 7    0.358  0.56

pair_modify shift yes

bond_coeff 1 400.0  1.00
special_bonds lj 1 1 1

angle_coeff 1 100.0 90.0

dihedral_coeff 1 100.0 -1 2
dihedral_coeff 2 100.0 -1 2


mass *  0.1667
mass 2  0.0812
mass 3  0.0812
mass 4  0.0812
mass 5  0.0812
mass 6  0.0812
mass 7  0.0812
mass 14 0.0008

run 0
velocity  all create $T  492459 rot yes mom yes units box
fix 1 all rigid/small molecule langevin $T $T 0.0001 428984
run 0
neigh_modify exclude molecule/intra all
comm_modify cutoff 3.5
velocity all scale $T
run 0

compute         c1  all cluster/atom 1.3
compute         maxcluid all reduce max c_c1
compute         cc1 all chunk/atom c_c1 compress yes
compute         size all property/chunk cc1 count

dump            1 all custom 50000 dump_${T}.lammpstrj id type x y z c_c1
dump_modify     1 append no sort id

timestep        1e-4
thermo          50000
thermo_style    custom step temp epair emol etotal press
run             5000000
write_data      last_conf_${T}
