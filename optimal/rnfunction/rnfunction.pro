######################################################################
# Automatically generated by qmake (3.0) Sat Sep 3 16:45:57 2016
######################################################################

TEMPLATE = lib
TARGET   = rnfunction
CONFIG  += console
CONFIG  -= app_bundle
CONFIG  -= qt
CONFIG  += dll
DESTDIR  = ../bin

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

DEFINES += MINIMUM_LIBRARY

# Input
HEADERS += bealesfunction.h boothfunction.h rosenbrock.h
SOURCES += bealesfunction.cpp boothfunction.cpp rosenbrock.cpp
