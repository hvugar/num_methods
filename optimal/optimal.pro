TEMPLATE = subdirs
CONFIG  += ordered
CONFIG  += -std=c++11
CONFIG  += c++11

SUBDIRS += minimum
#SUBDIRS += rnfunction
#SUBDIRS += control
#SUBDIRS += third-party
#SUBDIRS += point
#SUBDIRS += parabolic
#SUBDIRS += hyperbolic
#SUBDIRS += discrete
#SUBDIRS += imager
#SUBDIRS += imaging
SUBDIRS += border
#SUBDIRS += problem0H
#SUBDIRS += problem1
#SUBDIRS += problem1P
#SUBDIRS += problem2P
#SUBDIRS += problem1H
#SUBDIRS += problem3P
#SUBDIRS += problem2H
#SUBDIRS += heightMap3d
#SUBDIRS += kamil
SUBDIRS += main

problem1.depends = minimum
problem1P.depends = minimum
problem2P.depends = minimum
problem1H.depends = minimum
problem2H.depends = minimum
problem3P.depends = minimum
main.depends = minimum
main.depends = problem1
main.depends = problem1P
main.depends = problem2P
main.depends = problem1H
main.depends = problem2H
main.depends = problem3P
