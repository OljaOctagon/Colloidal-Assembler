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
	global iop
	
	
	graphics top delete all
	graphics top material AOShiny
	graphics top color orange
 	#graphics top color orange
	for {set icube 0} {$icube < $ncubes} {incr icube} {
		#set coords [[atomselect top "atomicnumber==$icube" frame $vmd_frame([molinfo top])] get {x y z}]
		#puts [lindex data $
		set istart [expr $icube*6]
		set iend [expr ($icube+1)*6]
		#set coords [[atomselect top "index>=$istart && index<$iend" frame $vmd_frame([molinfo top])] get {x y z}]
		#puts $coords



		set iline [expr $nlines*$vmd_frame([molinfo top])+$icube*6+2]
		set V0 [list [lindex $data [expr $iline + 0] 1] [lindex $data [expr $iline + 0] 2] [lindex $data [expr $iline + 0] 3]]
		set V1 [list [lindex $data [expr $iline + 1] 1] [lindex $data [expr $iline + 1] 2] [lindex $data [expr $iline + 1] 3]]
		set V2 [list [lindex $data [expr $iline + 2] 1] [lindex $data [expr $iline + 2] 2] [lindex $data [expr $iline + 2] 3]]
		set V3 [list [lindex $data [expr $iline + 3] 1] [lindex $data [expr $iline + 3] 2] [lindex $data [expr $iline + 3] 3]]
		set V4 [list [lindex $data [expr $iline + 4] 1] [lindex $data [expr $iline + 4] 2] [lindex $data [expr $iline + 4] 3]]
		set V5 [list [lindex $data [expr $iline + 5] 1] [lindex $data [expr $iline + 5] 2] [lindex $data [expr $iline + 5] 3]]
		
		set coords [list $V0 $V1 $V2 $V3 $V4 $V5]

		vmd_draw_cube top $coords
		
		
		set p1ax [ lindex $V0 0 ]
		set p1ay [ lindex $V0 1 ]
		set p1az [ lindex $V0 2 ]

		set p1bx [ lindex $V3 0 ]
		set p1by [ lindex $V3 1 ]
		set p1bz [ lindex $V3 2 ]

		set p2ax [ lindex $V1 0 ]
        set p2ay [ lindex $V1 1 ]
        set p2az [ lindex $V1 2 ]
                                          
		set p2bx [ lindex $V4 0 ]
    	set p2by [ lindex $V4 1 ]
		set p2bz [ lindex $V4 2 ] 

		set p3ax [ lindex $V2 0 ]  				
        set p3ay [ lindex $V2 1 ]
		set p3az [ lindex $V2 2 ]
                                           
        set p3bx [ lindex $V5 0 ]
        set p3by [ lindex $V5 1 ]
        set p3bz [ lindex $V5 2 ]
                                        
		set pp1 0.25

		set p1_x [ expr $p1ax + $pp1*($p1bx - $p1ax)]
		set p1_y [ expr $p1ay + $pp1*($p1by - $p1ay)]
		set p1_z [ expr $p1az + 0.5 *($p1bz - $p1az)]

		set p2_x [ expr $p1bx + $pp1*($p1ax - $p1bx)] 	
        set p2_y [ expr $p1by + $pp1*($p1ay - $p1by)]
        set p2_z [ expr $p1bz + 0.5 *($p1az - $p1bz)]

		set p3_x [ expr $p2ax + $pp1*($p2bx - $p2ax)]
        set p3_y [ expr $p2ay + $pp1*($p2by - $p2ay)]
        set p3_z [ expr $p2az + 0.5 *($p2bz - $p2az)]
                                                     
        set p4_x [ expr $p2ax + $pp1*($p2ax - $p2bx)] 
        set p4_y [ expr $p2ay + $pp1*($p2ay - $p2by)]
        set p4_z [ expr $p2az + 0.5 *($p2az - $p2bz)]

		set p5_x [ expr $p3ax + $pp1*($p3bx - $p3ax)]
        set p5_y [ expr $p3ay + $pp1*($p3by - $p3ay)]
        set p5_z [ expr $p3az + 0.5 *($p3bz - $p3az)]
                                                     
        set p6_x [ expr $p3ax + $pp1*($p3ax - $p3bx)] 
        set p6_y [ expr $p3ay + $pp1*($p3ay - $p3by)]
        set p6_z [ expr $p3az + 0.5 *($p3az - $p3bz)]

		graphics top color red
		set r_patch 0.05

		set p1_vec [list $p1_x $p1_y $p1_z ]
		set p2_vec [list $p2_x $p2_y $p2_z ]
		set p3_vec [list $p3_x $p3_y $p3_z ] 	    
	    set p4_vec [list $p4_x $p4_y $p4_z ]
	    set p5_vec [list $p5_x $p5_y $p5_z ]
	    set p6_vec [list $p6_x $p6_y $p6_z ]

	    graphics top color blue
		draw sphere $p1_vec radius $r_patch  resolution 30		
		draw sphere $p2_vec radius $r_patch  resolution 30		
		draw sphere $p3_vec radius $r_patch  resolution 30		
		draw sphere $p4_vec radius $r_patch  resolution 30		
		draw sphere $p5_vec radius $r_patch  resolution 30		
		draw sphere $p6_vec radius $r_patch  resolution 30		


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
set ncubes [expr [[atomselect top all frame $vmd_frame([molinfo top])] num]/6]

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

