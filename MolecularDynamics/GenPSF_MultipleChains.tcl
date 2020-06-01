################################################################################
# GenPSF multiple chains     #
##############################
# This script reads a .pdb file with 6 chains and saves each chain in a
# separate .pdb. Also generates .psf file for the whole protein and each chain.
#
# Written by Ana Paula Vargas R.
# Bioinformatics lab - UPCH. Lima, Peru.
# Last updated: 31-may-20.
###############################################################################


# Modelo 3PXG - All chains
#

# Cargar molécula
mol new modelo3pxg_all.pdb
set molID [molinfo top]

# Separar molécula por cadenas
set ChainA [atomselect $molID "chain A and protein"]
set ChainB [atomselect $molID "chain B and protein"]
set ChainC [atomselect $molID "chain C and protein"]
set ChainD [atomselect $molID "chain D and protein"]
set ChainE [atomselect $molID "chain E and protein"]
set ChainF [atomselect $molID "chain F and protein"]

# Guardar cadena en PDB distinto
$ChainA writepdb Chain_A.pdb
$ChainA delete
$ChainB writepdb Chain_B.pdb
$ChainB delete
$ChainC writepdb Chain_C.pdb
$ChainC delete
$ChainD writepdb Chain_D.pdb
$ChainD delete
$ChainE writepdb Chain_E.pdb
$ChainE delete
$ChainF writepdb Chain_F.pdb
$ChainF delete

mol delete $molID

# Procedimiento para crear archivos por cadena.
proc CreateChainFiles {segName PDBName OutName} {

# Generar PSF y PDB para cada cadena.

        package require psfgen; # Cargar psfgen. Borrar datos anteriores.
        resetpsf
        psfcontext reset

        topology 01_Topology/toppar/top_all36_prot.rtf; # Cargar topología.
        topology 01_Topology/toppar/toppar_all36_prot_modify_res.str

        pdbalias residue HIS HSE
	pdbalias atom ILE CD1 CD

# Crear segmento.
        segment $segName {
                pdb $PDBName
        }

        coordpdb $PDBName $segName
        guesscoord

# Escribir archivos.
        writepdb $OutName.pdb
        writepsf $OutName.psf

# Limpieza.
        resetpsf
        psfcontext reset
}

# Ejecutar procedimiento:
CreateChainFiles A Chain_A.pdb A
CreateChainFiles B Chain_B.pdb B
CreateChainFiles C Chain_C.pdb C
CreateChainFiles D Chain_D.pdb D
CreateChainFiles E Chain_E.pdb E
CreateChainFiles F Chain_F.pdb F

# Limpieza.
resetpsf
psfcontext reset



######################################
# Generar PSF para toda la molécula
######################################

package require psfgen
resetpsf

readpsf A.psf
coordpdb A.pdb

readpsf B.psf
coordpdb B.pdb

readpsf C.psf
coordpdb C.pdb

readpsf D.psf
coordpdb D.pdb

readpsf E.psf
coordpdb E.pdb

readpsf F.psf
coordpdb F.pdb


writepsf modelo3pxg_all_test1.psf
writepdb modelo3pxg_all_test1.pdb

resetpsf

exit
