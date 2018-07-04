TEMPLATE = subdirs
CONFIG += ordered

SUBDIRS += minimum \
    exporter
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

main.depends = minimum
