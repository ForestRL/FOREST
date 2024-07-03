# FOREST
FORcasting Events from Supernovae Theoretical modeling (FOREST)

### Ouput format
    Ev#, time[s], Ev ene[MeV], Nu ene[MeV], theta, phi, x, y, z, Ev id, FV

## Examples

### Example 1
    python3 src/forest.py -model analytic_rate -distance 10 -detector superk -pns_m 1.5 -pns_r 11.8 -etot 1e53 -gbeta 1.6  -end_time 200 -detector superk

### Example 2
    python3 src/forest.py -model numerical_spectra -spectra_file z9.6_ver2_90s.bary.1.36_grav.1.26.dat -distance 10  -detector superk
