# executables in $PREFIX/bin :

bin_PROGRAMS = gensolvent gmxtop2solute heatmap ts4sdump eprism3d
dist_bin_SCRIPTS = src/generate-idc.sh

gensolvent_SOURCES = src/gensolvent.cpp
gmxtop2solute_SOURCES = src/gmxtop2solute.cpp
heatmap_SOURCES = src/heatmap.cpp
ts4sdump_SOURCES = src/ts4sdump.cpp
eprism3d_SOURCES = src/eprism3d.cpp

LDADD = $(EXTRA_LDFLAGS) $(LDFLAGS)
eprism3d_CXXFLAGS = $(CFLAGS) $(EXTRA_CFLAGS)
gensolvent_CXXFLAGS = $(CFLAGS) $(EXTRA_CFLAGS)

# solvent and solute files in $PREFIX/share/eprism :

solventdir = ${pkgdatadir}/solvent
solvent_DATA = solvent/gvv.tip3p.1DDRISM.01A.txt solvent/gvv.tip3p.1DDRISM.025A.txt solvent/tip3p-amber14.01A.gaff solvent/tip3p-amber14.025A.gaff

solutedir = ${pkgdatadir}/solute
solute_DATA = solute/methane.gro solute/methane.pdb solute/methane.prmtop solute/methane.top solute/methane.solute solute/methane.xtc

EXTRA_DIST = ${solvent_DATA} ${solute_DATA}

