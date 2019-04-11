# This tcl script loads a trajectory of liquid water simulated
# with NAMD at 300K into VMD, and then calls the script pdf.tcl to 
# calculate the pair distribution function.
#
# Launch VMD, open the Tk Console and then type "source calcpdf.tcl"
# or
# type "vmd -dispdev text -e calcpdf.tcl" in the unix command line

# Load the psf file into VMD
mol load psf ICES.psf

# Load the trajectory
animate read dcd melt-300-01.dcd skip 10 waitfor all

# Select oxygen atoms of water 
set ox [atomselect top "name OH2"]

# Calls pdf.tcl
# Parameters for pdf: 
# molecule id, selection, Size of Box along X axis, along Y axis, along Z 
# axis, Number of bins, size of bins (note that Number of bins * size of 
# bins gives you the range for which the pair distribution function is 
# computed. 
source pdf.tcl
pdf top $ox 23.623 22.406 27.1759 100 0.1

