################################################################################
# Photos RNA            #
#########################
# This script loads a molecule in VMD, rotates it and takes a number of pictures
# ready for publication.
# This script works, but needs to be setted manually and
# is currently being re-written for full automatization.
#
# Written by Ana Paula Vargas.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 17-feb-20.
###############################################################################




###########################################################
# PROCEDURE: Take X photos of protein in different angles #
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
# actualizado 03/dic/19
# Copiar y pegar los dos proc y usar photos_rna pnum Deg
# por ejemplo: photos_rna 8 45 tomará 8 fotos, una cada 45 grados
# Solo está implementada la rotación en el eje Y
# Falta modificar para que el programa use una RES y VEC indicada por el usuario
# Solo implementado en Tachyon




######## ROTATE AXIS #########
# Proc. by Andrew Dalke
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


######## RENDER #########

proc photos_rna {pnum Deg} {
	set c 0
	while {$c < $pnum} {
		set outName "Fig6_$c"
		render Tachyon $outName "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %s -format TARGA -res 8000 6400 -o %s.tga
		rotate_axis {0 1 0} $Deg
		set c [expr {$c + 1}]
	}
}



##########################################
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
