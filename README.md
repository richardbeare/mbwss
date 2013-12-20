mbwss
=====

Marker based watershed scalper

This repository contains code for two tools for extracting brains in
T1 MRI scans - one for humans and one for macaques. The methods are
being published in Frontiers in Neuroinformatics, special issue on
ITK.

Early versions of repository are for the Oxford guys to try out the
macaque brain extraction tool.

Building - the easy way
--------

Prerequisites : cmake 2.8.7, git

git clone https://github.com/richardbeare/mbwss.git

cd mbwss

git submodule init

git submodule update

cd ..

mkdir BuildMBWSS

cd BuildMBWSS

ccmake ../mbwss/SuperBuild

make

An executable named scalperWSMacaque will be generated in 

MBWSS-build/bin/scalperWSMacaque


The steps above will download and install ITK 4.4. If you already
have it installed, do this instead:

ccmake ../mbwss/

Set the location of ITK4.4

make

Citations
---------

If you find these tools useful in your research, please cite:

"Brain extraction using the watershed transform from markers., R. Beare, J. Chen, C.L. Adamson, T. Silk, D.K. Thompson, J.Y.M. Yang, V.A. Anderson, M.L. Seal and A.G. Wood. Frontiers in Neuroinformatics, 7(32), December 2013."


@Article{MBWSSBeare,
  author =       {Beare, R. and Chen, J. and Adamson, C.L. and  Silk, T. and Thompson, D.K. and Yang, J.Y.M. and Anderson, V.A. and Seal, M.L. and Wood, A.G.},
  title =        {Brain extraction using the watershed transform from markers.},
  journal =      {Frontiers in Neuroinformatics},
  year =         2013,
  volume =       7,
  number =       32,
  month =        {December},
  doi = {doi:10.3389/fninf.2013.00032}
}


Notes on the macaque brain extractor:
-------------------------------------

Name: scalperWSMacaque

This tool isn't yet well tested. The main issue with the scans I've
developed for is extremely heavy brightness inhomogeneity. There is a
built in (optional) bias correction step that is intended to be basic
and fast and normalise the intensity sufficiently for edges to be
detectable. It is built around largish mean filters, the size of which
is the only exposed parameter. I've had good results so far setting
the size to 10mm, which makes the very faint areas visible again. If
results are very bad, this is probably a productive place to start.

Scans need to be in approximately MNI orientation - i.e RL PA IS (or
LR PA IS). Angle isn't particulary critical.

Usage:

scalperWSMacaque --help

The call I've been using for the scans I've tested so far is:

scalperWSMacaque -i inputfile.nii.gz -o /tmp/output --refine --biascorrect --preopening 1 --smoothGradSigma 1.0

I've found a couple of common sources of error. One is strong inhomogeneity 
near image edges and the other is some of the optic nerves get partially
segmented. The latter problem can probably be reduced by increasing
the size of the smoothing, and it doesn't look very serious to me anyway.
The internal bias field corrector is very basic. I've changed it recently
to match the one I'm using for human scans. If there are
problems with new scans then perhaps we need to switch back to the original version.

