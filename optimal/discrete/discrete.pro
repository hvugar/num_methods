######################################################################
# Automatically generated by qmake (3.0) Sat Sep 3 16:41:26 2016
######################################################################

TEMPLATE = lib
TARGET   = discrete
CONFIG  += console
CONFIG  -= app_bundle
CONFIG  -= qt
CONFIG  += dll
DESTDIR  = ../bin

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

# Input
HEADERS += discreteheat.h discretehyperbolic.h discretehyperbolic1.h
SOURCES += discreteheat.cpp discretehyperbolic.cpp discretehyperbolic1.cpp
