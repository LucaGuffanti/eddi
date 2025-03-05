# Change the paths to the data files here
# The paths will be used to load the data in the AtomDescriptor class

from utils.download_data import load_vdw_radii

PERIODIC_TABLE_PATH = 'data/periodic_table.csv'
VDW_RADII_PATH = 'data/vdw_radii.csv'
data = load_vdw_radii(vdw_radii_path=VDW_RADII_PATH, periodic_table_path=PERIODIC_TABLE_PATH, verbose=True)

