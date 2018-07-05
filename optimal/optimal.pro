TEMPLATE = subdirs
CONFIG += ordered

SUBDIRS += minimum
#SUBDIRS += rnfunction
#SUBDIRS += control
#SUBDIRS += third-party
#SUBDIRS += point
#SUBDIRS += border
#SUBDIRS += parabolic
#SUBDIRS += hyperbolic
#SUBDIRS += discrete
#SUBDIRS += imaging
#SUBDIRS += imager
SUBDIRS += problem2H
SUBDIRS += main
SUBDIRS += exporter

main.depends = minimum
problem2H.depends = minimum
exporter.depends = problem2H
