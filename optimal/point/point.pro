######################################################################
# Automatically generated by qmake (3.0) Sat Sep 3 16:49:25 2016
######################################################################

TEMPLATE = lib
TARGET   = point
CONFIG  += console
CONFIG  -= app_bundle
CONFIG  -= qt
CONFIG  += dll
DESTDIR  = ../bin

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum


# Input
HEADERS += pointcontrol.h pointcontrol1.h pointcontrol11.h pointcontrol2.h
SOURCES += pointcontrol.cpp \
           pointcontrol1.cpp \
           pointcontrol11.cpp \
           pointcontrol2.cpp
