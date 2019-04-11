# Script to compute pair distribution function from a selection
# (Based on an example from "Computer Simulation of Liquids", 
# M. P. Allen and D. J. Tildesley)
#
# Example:
#
# Load a trajectory in VMD
# Then, in the Tk Console type: 
# source pdf.tcl
# pdf <mol> <sel> <BOXXL> <BOXYL> <BOXZL> <maxbin> <delr> 
#
# <mol>: Molecule id (0, 1, 2... or top)
# <sel>: A selection (created with the command atomselect)
# <BOXXL>: Size of the periodic box along the X axis
# <BOXYL>: Size of the periodic box along the Y axis
# <BOXZL>: Size of the periodic box along the Z axis
# <maxbin>: Number of bins used to compute the pair distribution function
# <delr>: Size of the bins
# Note that <maxbin> x <delr> = range over which the pair distribution 
# function is computed.
#
# Marcos Sotomayor
# sotomayo@ks.uiuc.edu

# For a given selection, computes the pair distribution function
# returns a list with the computed values
proc pdf {mol sel BOXXL BOXYL BOXZL maxbin delr} {

   
    for {set k 0} {$k < $maxbin} {incr k} {
	lappend hist 0
    }
 
    set ind [list [$sel get index]]
    set indl [llength [lindex $ind 0]]

    set NF [molinfo $mol get numframes]
    for {set m 0} {$m < [expr $NF]} {incr m} {
	puts "frame: $m/$NF"
	$sel frame $m
	for {set i 0} {$i < [expr $indl-1]} {incr i} {
	    for {set j [expr $i+1]} {$j < [expr $indl]} {incr j} {

		set atom1 [atomselect $mol "index [lindex $ind 0 $i]"]
		set atom2 [atomselect $mol "index [lindex $ind 0 $j]"]
	
		set rxij [expr [$atom1 get x] - [$atom2 get x]]
		set rxijm [expr $rxij - $BOXXL*round(double($rxij)/$BOXXL)]
		set ryij [expr [$atom1 get y] - [$atom2 get y]]
		set ryijm [expr $ryij - $BOXYL*round(double($ryij)/$BOXYL)]
		set rzij [expr [$atom1 get z] - [$atom2 get z]]
		set rzijm [expr $rzij - $BOXZL*round(double($rzij)/$BOXZL)]
		set rij [expr sqrt(pow($rxijm,2)+pow($ryijm,2)+pow($rzijm,2))]
		set bin [expr round(double($rij)/$delr)]
		if {$bin < $maxbin} then {
		    lset hist $bin [expr [lindex $hist $bin]+2]
		}
		$atom1 delete
		$atom2 delete
	    } 
	}
    }

    #for {set i 0} {$i< $maxbin} {incr i} {
	#puts "[lindex $hist $i]"
    #} 

    set rho [expr double([$sel num])/($BOXXL*$BOXYL*$BOXZL)]
    for {set i 0} {$i< $maxbin} {incr i} {
	set rlower [expr ($i-1)*$delr]
        set rupper [expr $rlower + $delr]
	set nideal [expr 4.0*3.141592654*$rho/3.0*(pow($rupper,3)-pow($rlower,3))]
	puts "[expr $i*$delr] [expr double([lindex $hist $i])/$NF/[$sel num]/$nideal]"
    }
    return "done"
}

# Do average and normalization
