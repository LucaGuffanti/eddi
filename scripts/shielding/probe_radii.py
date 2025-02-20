# probe_radii.py
#
# This scripts approximates the atomic radius of an atom given its atomic number, using the following approach:
# 1. The radial component of the electronic wavefunction is approximated via Slater's rules.
# 2. The electron density is computed around the atom by summing the modules squared of the wavefunctions for all the electrons.
# 3. As the integral in space of the electron density is equal to the number of electrons, the radius is computed as the radius of a sphere that
#    encapsulates, up to a user-chosen threshold, the number of electrons in the atom.