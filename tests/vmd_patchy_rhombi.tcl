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
		
		set p1ax [ lindex $V0 0 ]
		set p1ay [ lindex $V0 1 ]
		set p1az [ lindex $V0 2 ]

		set p1bx [ lindex $V7 0 ]
		set p1by [ lindex $V7 1 ]
		set p1bz [ lindex $V7 2 ]

		set p2ax [ lindex $V2 0 ]
                set p2ay [ lindex $V2 1 ]
                set p2az [ lindex $V2 2 ]
                                          
	        set p2bx [ lindex $V7 0 ]
                set p2by [ lindex $V7 1 ]
		set p2bz [ lindex $V7 2 ] 

		set p3ax [ lindex $V0 0 ]  				
                set p3ay [ lindex $V0 1 ]
                set p3az [ lindex $V0 2 ]
                                           
                set p3bx [ lindex $V5 0 ]
                set p3by [ lindex $V5 1 ]
                set p3bz [ lindex $V5 2 ]
                                           
                set p4ax [ lindex $V1 0 ]
                set p4ay [ lindex $V1 1 ]
                set p4az [ lindex $V1 2 ]
                                          
                set p4bx [ lindex $V6 0 ]
                set p4by [ lindex $V6 1 ]
                set p4bz [ lindex $V6 2 ] 


		set pp1 0.2
		set pp2 0.8

		set p1_x [ expr $p1ax + $pp1*($p1bx - $p1ax)]
		set p1_y [ expr $p1ay + $pp1*($p1by - $p1ay)]
		set p1_z [ expr $p1az + 0.5*($p1bz - $p1az)]

		set p2_x [ expr $p2ax + $pp2*($p2bx - $p2ax)] 	
        set p2_y [ expr $p2ay + $pp2*($p2by - $p2ay)]
        set p2_z [ expr $p2az + 0.5*($p2bz - $p2az)]

		set p3_x [ expr $p3ax + $pp2*($p3bx - $p3ax)]
        set p3_y [ expr $p3ay + $pp2*($p3by - $p3ay)]
        set p3_z [ expr $p3az + 0.5*($p3bz - $p3az)]
                                                     
        set p4_x [ expr $p4ax + $pp2*($p4bx - $p4ax)] 
        set p4_y [ expr $p4ay + $pp2*($p4by - $p4ay)]
        set p4_z [ expr $p4az + 0.5*($p4bz - $p4az)]

	
		graphics top color red
		set r_patch 0.1

		set p1_vec [list $p1_x $p1_y $p1_z ]
		set p2_vec [list $p2_x $p2_y $p2_z ]
		set p3_vec [list $p3_x $p3_y $p3_z ] 	    
	    set p4_vec [list $p4_x $p4_y $p4_z ]

	    graphics top color blue
		draw sphere $p1_vec radius $r_patch  resolution 30		
		draw sphere $p2_vec radius $r_patch  resolution 30		
		graphics top color red 
		draw sphere $p3_vec radius $r_patch  resolution 30	
        draw sphere $p4_vec radius $r_patch  resolution 30	
	
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

