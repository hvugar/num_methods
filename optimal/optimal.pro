TEMPLATE = subdirs
CONFIG += ordered

SUBDIRS += minimum
#SUBDIRS += rnfunction
#SUBDIRS += control
#SUBDIRS += third-party
#SUBDIRS += point
#SUBDIRS += parabolic
#SUBDIRS += hyperbolic
#SUBDIRS += discrete
#SUBDIRS += imager
SUBDIRS += imaging
#SUBDIRS += border
#SUBDIRS += problem1
#SUBDIRS += problem2P
SUBDIRS += problem2H
SUBDIRS += heightMap3d
SUBDIRS += main

main.depends = minimum
main.depends = problem1
main.depends = problem2P
main.depends = problem2H
problem1.depends = minimum
problem2P.depends = minimum
problem2H.depends = minimum
