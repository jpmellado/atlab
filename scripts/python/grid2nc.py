import sys
import atlab

# getting data from stdin
if len(sys.argv) <= 1:
    print("Usage: python $0 list-of-grid-files.")
    quit()

files = sorted(sys.argv[1:])

###########################################################
# process grid files
for file in files:
    print("Processing file %s ..." % file)

    grid = atlab.Grid()
    grid.read(filename=file)

    grid.writeNetCDF(filename=file)