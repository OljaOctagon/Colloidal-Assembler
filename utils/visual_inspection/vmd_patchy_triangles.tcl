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
	set c010 [lindex $coords 2]
	set c001 [lindex $coords 3]
	set c101 [lindex $coords 4]
	set c011 [lindex $coords 5]


	# top and bottom surface
	graphics top triangle $c000 $c100 $c010
	graphics top triangle $c001 $c101 $c011

	# sides bottom left triangle
	graphics top triangle $c000 $c100 $c001
	graphics top triangle $c100 $c010 $c101
	graphics top triangle $c010 $c000 $c011

	# sides upper right triangle
	graphics top triangle $c001 $c100 $c101
	graphics top triangle $c101 $c010 $c011
	graphics top triangle $c011 $c000 $c001

}

# DELETE OLD CUBES AND DRAW NEW CUBES IN CURRENT FRAME
proc drawpolyframe { name element op } {

	global vmd_frame
	global ncubes
	global nlines
	global data
	global data_op1
	global data_op2
	global iop
	global nlines_op
	
	
	graphics top delete all
	graphics top material AOShiny
	graphics top color orange
 	#graphics top color orange
	for {set icube 0} {$icube < $ncubes} {incr icube} {
		#set coords [[atomselect top "atomicnumber==$icube" frame $vmd_frame([molinfo top])] get {x y z}]
		#puts [lindex data $
		set istart [expr $icube*8]
		set iend [expr ($icube+1)*8]
		#set coords [[atomselect top "index>=$istart && index<$iend" frame $vmd_frame([molinfo top])] get {x y z}]
		#puts $coords
		set iline_op [expr $nlines_op*$vmd_frame([molinfo top])+$icube+2]
		set op1 [lindex $data_op1 $iline_op $iop]

		if {($op1 == 1 ) } {
			graphics top color 30
		}

		if {($op1 == 2 ) } {
			graphics top color blue2
		}

		if {($op1 == 3) } {
	       		graphics top color 31
		}

		set iline [expr $nlines*$vmd_frame([molinfo top])+$icube*6+2]
		set V0 [list [lindex $data [expr $iline + 0] 1] [lindex $data [expr $iline + 0] 2] [lindex $data [expr $iline + 0] 3]]
		set V1 [list [lindex $data [expr $iline + 1] 1] [lindex $data [expr $iline + 1] 2] [lindex $data [expr $iline + 1] 3]]
		set V2 [list [lindex $data [expr $iline + 2] 1] [lindex $data [expr $iline + 2] 2] [lindex $data [expr $iline + 2] 3]]
		set V3 [list [lindex $data [expr $iline + 3] 1] [lindex $data [expr $iline + 3] 2] [lindex $data [expr $iline + 3] 3]]
		set V4 [list [lindex $data [expr $iline + 4] 1] [lindex $data [expr $iline + 4] 2] [lindex $data [expr $iline + 4] 3]]
		set V5 [list [lindex $data [expr $iline + 5] 1] [lindex $data [expr $iline + 5] 2] [lindex $data [expr $iline + 5] 3]]
		
		set coords [list $V0 $V1 $V2 $V3 $V4 $V5]

		vmd_draw_cube top $coords
		
	
		graphics top color orange 
		#graphics top color orange

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


# LOAD USER DATA  FOR COLORING OP1
set file "color_op.dat"
set fp [open $file r]
set data_op1 [split [read $fp] "\n"]; puts ""
close $fp
unset fp
set nlines_op [expr [lindex $data_op1 0 0]+2]


# GOTO FIRST FRAME
set vmd_frame([molinfo top]) 0
color Display Background white 
display projection orthographic
display depthcue off
display resize 2000 1800 
display rendermode Acrobat3D
axes location off

