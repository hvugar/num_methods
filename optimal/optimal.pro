TEMPLATE = subdirs

SUBDIRS += minimum \
    third-party
#SUBDIRS += imaging
#SUBDIRS += imager
SUBDIRS += rnfunction
#SUBDIRS += control
#SUBDIRS += point
SUBDIRS += border
SUBDIRS += parabolic
#SUBDIRS += hyperbolic
#SUBDIRS += discrete
SUBDIRS += main

CONFIG += ordered
