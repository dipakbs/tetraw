set pdb_file [lindex $argv 0]
set num_atoms [lindex $argv 1]
mol new $pdb_file 
set oxy [atomselect top "name O"]
set index_o [$oxy get index]

set file [open input_${num_atoms}.dat w]
foreach o $index_o {
set water [atomselect top "index $o"]
set coord_x [$water get {x}]
set coord_y [$water get {y}]
set coord_z [$water get {z}]
puts $file "$o,$coord_x,$coord_y,$coord_z"
$water delete
}
close $file
exit
