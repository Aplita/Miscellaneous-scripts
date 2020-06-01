################################################################################
# Photos RNA - WIP      #
#########################
# This script loads a molecule in VMD, rotates it and takes a number of pictures
# ready for publication.
# This script works, but needs to be setted manually and
# is currently being re-written for full automatization.
# THIS IS ITS LAST VERSION (not tested or finished yet)
#
# Written by Ana Paula Vargas.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 17-feb-20.
###############################################################################





###########################################################
# SCRIPT: Take X photos of protein in different angles #
###########################################################

# vars to set:
# vec: a vector of shape {x y z} filled with 1/0
# to toggle rotation in corresponding axis.
# deg: degree of rotation.
# pnum: number of photos to take.
# c: internal counter.
# resX: resolution X axis.
# resY: resolution Y axis.
#

# última versión:
# actualizado 31/ene/20
# Copiar y pegar los dos proc y usar photos_rna pnum Deg
# por ejemplo: photos_rna 8 45 tomará 8 fotos, una cada 45 grados
# Solo está implementada la rotación en el eje Y
# Falta modificar para que el programa use una RES y VEC indicada por el usuario
# Solo implementado en Tachyon


puts "This script will rotate your molecule and render (via Tachyon) a HD picture in each pose.\n
Last updated: 17 feb 20\n
Author: https://github.com/Aplita\n_______________________________"
puts ""

puts "Defaults:\n
1) Resolution: 8000x6400\n
2) Number of photos: 4\n
3) Rotation: only Y axis\n
4) Rotation degree: 90°\n
5) Output Name: Fig_X \n
(X changes according to pic number, as in Fig_1, Fig_2, etc.)
\n_______________________________"
puts ""


set C 0
set c 0

while {$C == 0} {
	puts -nonewline "Keep defaults? (y/n): "
	gets stdin answ
	if {$answ == "Y" | $answ == "y"} {
		set C 1
	} elseif {$answ == "N" | $answ == "n"} {
		set C 2
	} else {
		puts "Please answer as y (yes) or n (no)"
	}
}


###########################################################
# PROCEDURES						  #
###########################################################


######## ROTATE AXIS ##########
## procedure by Andrew Dalke ##

proc rotate_axis {vec deg {molid top}} {
    # get the current matrix
    lassign [molinfo $molid get rotate_matrix] curr
    # the transformation matrix
    set r [trans axis $vec $deg]
    # get the new matrix
    set m [transmult $r $curr]
    # and apply it to the molecule
    molinfo $molid set rotate_matrix "{ $m }"
}


######## RENDER AND ROTATE #########
proc photos_rna {Pnum Deg Res Vec OutName} {
        set c 0
        while {$c < $Pnum} {
                render Tachyon $OutName "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -res $Res -o %s.tga
                rotate_axis $Vec $Deg
                set c [expr {$c + 1}]
        }
}

####### INPUT VERIFICATION ########
proc input_ver {Inp} {
	set L [llength [split $Inp x]]
	if {$L == 2} {
		return "T"
	} else {
		return "F"
	}
}

###### SET ROTATION DEGREE ########
proc set_rot_deg {Vec} {
	set tX [lindex [split $Vec " "] 0]
        set tY [lindex [split $Vec " "] 1]
        set tY [lindex [split $Vec " "] 2]

	if {$tX == 1} {
		if {$tY == 1} {
			if {$tZ == 1} {
				# XYZ plane
				puts -nonewline "Degree of rotation in the XYZ plane: "
				gets stdin ndeg
			} else {
				# XY plane
				puts -nonewline "Degree of rotation in the XY plane: "
				gets stdin ndeg
			}
		} elseif {$tZ == 1} {
			# XZ plane
			puts -nonewline "Degree of rotation in the XZ plane: "
			gets stdin ndeg
		} else {
			# Only X axis
			puts -nonewline "Degree of rotation in X axis: "
			gets stdin ndeg
		}
	} elseif {$tY == 1} {
		if {$tZ == 1} {
			# YZ plane
			puts -nonewline "Degree of rotation in YZ plane: "
			gets stdin ndeg
		} else {
			# Only Y axis
			puts -nonewline "Degree of rotation in Y axis: "
			gets stdin ndeg
		}
	} else {
		# Only Z axis
		puts -nonewline "Degree of rotation in Z axis: "
		gets stdin ndeg
	}

	return $ndeg
}




###### SET ROTATION AXIS ##########
proc set_rot_ax {} {
	puts "Choose from the following options:"
	puts "1) Only X axis."
	puts "2) Only Y axis (default)."
	puts "3) Only Z axis."
	puts "4) X and Y."
	puts "5) X and Z."
	puts "6) Y and Z."
	puts "7) All (X, Y, Z)."
	puts -nonewline "Number of choice: "
	gets stdin inp

	if {$inp == 1} {
		set vect {1 0 0}
		set newD 90
	} elseif {$inp == 2} {
		set vect {0 1 0}
		set newD 90
	} elseif {$inp == 3} {
		set vect {0 0 1}
		set newD 90
	} elseif {$inp == 4} {
		# Need deg
		set vect {1 1 0}
		set newD [set_rot_deg $vect]
	} elseif {$inp == 5} {
		set vect {1 0 1}
		set newD [set_rot_deg $vect]
	} elseif {$inp == 6} {
		set vect {0 1 1}
		set newD [set_rot_deg $vect]
	} elseif {$inp == 7} {
		set vect {1 1 1}
		set newD [set_rot_deg $vect]
	} else {
		set newD 90
		puts "Bad input or out of range. Set to default."
	}
	return [list $vect $newD]
}


###########################################################
# EXECUTION						  #
###########################################################

if {$C == 1} {
	# Default options
	set deg 90
	set pnum 4
	set res "8000 6400"
	set outName "Fig_$c"
	set vec "{0 1 0}"
	puts "Working with default options..."
	photos_rna $pnum $deg $res $vec $outName
} else {
	set C 0
	while {C == 0} {
		puts -nonewline "Select number of option you would like to change: "
		gets stdin answ
		if {$answ == 1} {
			# Resolution
			puts -nonewline "Set new resolution (input example: 123x456): "
			gets stdin inp
			set nval [input_ver $inp]
			if {$nval == "T"} {
				set res [llength [split $inp x]]
				set C 1
			} else {
				puts "Invalid input. Input example: 123x456.\nExiting program."
				break
			}
		} elseif {$answ == 2} {
			# Number of pictures
			puts -nonewline "Set new number of pictures: "
			gets stdin inp
			set nval [string is entier $inp]
			if {$nval == 1} {
				# Is integer
				set pnum $answ
				set C 1
			} else {
				puts "Invalid input. Exiting program."
				break
			}
		} elseif {$answ == 3} {
			# Rotation axis
			# If user chooses more than 1 axis, program will ask for degree
			set inp [set_rot_ax]
			set vec [lindex $inp 0]
			set deg [lindex $inp 1]
			set C 1
		} elseif {$answ == 4} {
			# Rotation degree
		} elseif {$answ == 5} {
			# Output Name
		} else {
			puts "Invalid input. Please select a number according to the option you would like to change."
		}
	}

}
