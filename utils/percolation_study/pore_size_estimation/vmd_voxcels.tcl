# FUNCTION TO ENABLES TRACING
proc enabletrace {} {

	global vmd_frame
	trace variable vmd_frame([molinfo top]) w drawpolyframe

}

# FUNCTION TO DISABLE TRACING
proc disabletrace {} {

	global vmd_frame
	trace vdelete vmd_frame([molinfo top]) w drawpolyframe

}

# DRAW ONE CUBE
proc vmd_draw_cube {mol coords} {
	
	set c000 [lindex $coords 0]
	set c100 [lindex $coords 1]
	set c110 [lindex $coords 2]
	set c010 [lindex $coords 3]
	set c001 [lindex $coords 4]
	set c101 [lindex $coords 5]
	set c111 [lindex $coords 6]
	set c011 [lindex $coords 7]

	# top and bottom surface
	graphics top triangle $c000 $c100 $c110
	graphics top triangle $c000 $c010 $c110
	graphics top triangle $c001 $c101 $c111
	graphics top triangle $c001 $c011 $c111

	# left and right surface
	graphics top triangle $c000 $c100 $c101
	graphics top triangle $c000 $c001 $c101
	graphics top triangle $c010 $c110 $c111
	graphics top triangle $c010 $c011 $c111

	# front and back surface
	graphics top triangle $c000 $c010 $c011
	graphics top triangle $c000 $c001 $c011
	graphics top triangle $c100 $c110 $c111
	graphics top triangle $c100 $c101 $c111

}

# DELETE OLD CUBES AND DRAW NEW CUBES IN CURRENT FRAME
proc drawpolyframe { name element op } {

	global vmd_frame
	global ncubes
	global nlines
	global data
	global iop
	
	graphics top delete all
	graphics top material AOShiny
	graphics top color orange
	for {set icube 0} {$icube < $ncubes} {incr icube} {
		set istart [expr $icube*8]
		set iend [expr ($icube+1)*8]

		set iline [expr $nlines*$vmd_frame([molinfo top])+$icube*8+2]
		set V0 [list [lindex $data [expr $iline + 0] 1] [lindex $data [expr $iline + 0] 2] [lindex $data [expr $iline + 0] 3]]
		set V1 [list [lindex $data [expr $iline + 1] 1] [lindex $data [expr $iline + 1] 2] [lindex $data [expr $iline + 1] 3]]
		set V2 [list [lindex $data [expr $iline + 2] 1] [lindex $data [expr $iline + 2] 2] [lindex $data [expr $iline + 2] 3]]
		set V3 [list [lindex $data [expr $iline + 3] 1] [lindex $data [expr $iline + 3] 2] [lindex $data [expr $iline + 3] 3]]
		set V4 [list [lindex $data [expr $iline + 4] 1] [lindex $data [expr $iline + 4] 2] [lindex $data [expr $iline + 4] 3]]
		set V5 [list [lindex $data [expr $iline + 5] 1] [lindex $data [expr $iline + 5] 2] [lindex $data [expr $iline + 5] 3]]
		set V6 [list [lindex $data [expr $iline + 6] 1] [lindex $data [expr $iline + 6] 2] [lindex $data [expr $iline + 6] 3]]
		set V7 [list [lindex $data [expr $iline + 7] 1] [lindex $data [expr $iline + 7] 2] [lindex $data [expr $iline + 7] 3]]
		
		set coords [list $V0 $V1 $V2 $V3 $V4 $V5 $V6 $V7]

		vmd_draw_cube top $coords
		graphics top color orange 

		unset istart
		unset iline
		unset V0
		unset V1
		unset V2
		unset V3
		unset V4
		unset V5
		unset V6
		unset V7
		unset coords
	}

}

# GET FILE NAME FROM COMMAND LINE
set fname [lindex $argv 0]

# READ IN COORDINATES FROM FILE
mol new $fname type xyz waitfor all step 1

# DELETE AUTOMATICALLY GENERATED REPRESENTATION
mol delrep 0 top 

# GET NUMBER OF CUBES IN LAST FRAME
set ncubes [expr [[atomselect top all frame $vmd_frame([molinfo top])] num]/8]

# ENABLE TRACING
enabletrace

set iop 0

set fp [open $fname r]
set data [split [read $fp] "\n"]; puts ""
close $fp
unset fp
set nlines [expr [lindex $data 0 0]+2]

# GOTO FIRST FRAME
set vmd_frame([molinfo top]) 0
color Display Background white 
display projection orthographic
display depthcue off
display resize 2000 1800 
display rendermode Acrobat3D
axes location off

