set all [atomselect top all]

# determine the dimensions of the water box
set mes [measure minmax $all]
set xmin [lindex [lindex $mes 0] 0]
set ymin [lindex [lindex $mes 0] 1]
set zmin [lindex [lindex $mes 0] 2]
set xmax [lindex [lindex $mes 1] 0]
set ymax [lindex [lindex $mes 1] 1]
set zmax [lindex [lindex $mes 1] 2]

set O [atomselect top "name OH2"]

# set the number of atoms for which you want to count H-bonds
set atomcount 300.0
set totalbondnum 0
for {set i 1} {$i <= $atomcount} {incr i} {
  # pick a random water molecule (oxygen atom) and record its position
  set rand [expr int(rand()*[expr [$O num] - 1])]
  set oatom [atomselect top "index [lindex [$O get index] $rand]"]
  set xpos [lindex [lindex [$oatom get {x y z}] 0] 0]
  set ypos [lindex [lindex [$oatom get {x y z}] 0] 1]
  set zpos [lindex [lindex [$oatom get {x y z}] 0] 2]

  set buffer 3
  # the selected atom should be away from the surface (by <buffer> Angstroms)
  if {$xpos < [expr $xmin+$buffer] || $xpos > [expr $xmax-$buffer] || $ypos < [expr $ymin+$buffer] || $ypos > [expr $ymax-$buffer] || $zpos < [expr $zmin+$buffer] || $zpos > [expr $zmax-$buffer]} then {
    set i [expr $i-1]
  } else {
    # set the distance of a H-bond to 3.0 Angstroms and the angle to 30 degrees (default values)
    set dist 3.0
    set ang 30.0

    set temp1 [atomselect top "same residue as index [$oatom get index]"]
    set temp2 [atomselect top "same residue as exwithin 10 of index [$temp1 get index]"]
    # count the number of H-bonds for th molecule

    set count1 [llength [lindex [measure hbonds $dist $ang $temp1 $temp2] 0]]
    set count2 [llength [lindex [measure hbonds $dist $ang $temp2 $temp1] 0]]

    # and add it to the total so far
    set totalbondnum [expr $totalbondnum+$count1 + $count2]
    }
  }

puts "Each water molecule makes [expr double($totalbondnum/$atomcount)] hydrogen bonds on average with its neighbors."
