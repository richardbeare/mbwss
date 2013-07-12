mbwss
=====

Marker based watershed scalper

This repository contains code for two tools for extracting brains in
T1 MRI scans - one for humans and one for macaques. The methods are
being published in Frontiers in Neuroinformatics, special issue on
ITK.

Early versions of repository are for the Oxford guys to try out the
macaque brain extraction tool.

Building
--------

Prerequisites : cmake 2.8.7, git

git clone https://github.com/richardbeare/mbwss.git

mkdir BuildMBWSS
cd BuildMBWSS

ccmake ../mbwss/SuperBuild
make

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
results are very bad, this is probably a producting place to start.

Scans need to be in approximately MNI orientation - i.e RL PA IS (or
LR PA IS). Angle isn't particulary critical.

Usage:

scalperWSMacaque --help

