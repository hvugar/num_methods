TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH = ../minimum
LIBS += -L../bin -lminimum

INCLUDEPATH = ../cminimum
#LIBS += -L../bin -lcminimum

DESTDIR += ../bin

#DEFINES += MINIMUM_LIBRARY
#include(../minimum/minimum.pri)

SOURCES += main.cpp \
    utils.cpp \
    rosenbrock.cpp \
    bealesfunction.cpp \
    boothfunction.cpp \
    cfunction.cpp \
    cfunction1.cpp \
    cfunction2.cpp \
    cfunction3.cpp \
    pointcontrol.cpp \
    pointcontrol1.cpp \
    pointcontrol2.cpp \
    parabolic/1d/heatcontrol.cpp \
    parabolic/1d/heatcontroldeltaf.cpp \
    parabolic/1d/heatcontroldeltax.cpp \
    parabolic/2d/heatcontrol2d.cpp \
    parabolic/2d/heatcontrol2delta.cpp \
    parabolic/2d/heatcontrol2deltaf.cpp \
    parabolic/2d/heatcontrol2deltax.cpp \
    parabolic/1d/heatcontrol1.cpp \
    hyperbolic/hyperbolic2dx.cpp \
    hyperbolic/hyperbolic1dx.cpp \
    hyperbolic/hyperboliccontrol1d.cpp \
    hyperbolic/hyperboliccontrol1d2.cpp \
    hyperbolic/hyperboliccontrol1d3.cpp \
    hyperbolic/hyperboliccontrol1dt.cpp \
    hyperbolic/hyperboliccontrol1d4.cpp

HEADERS += \
    utils.h \
    rosenbrock.h \
    bealesfunction.h \
    boothfunction.h \
    cfunction.h \
    cfunction1.h \
    cfunction2.h \
    cfunction3.h \
    pointcontrol.h \
    pointcontrol1.h \
    pointcontrol2.h \
    parabolic/1d/heatcontrol.h \
    parabolic/1d/heatcontroldeltaf.h \
    parabolic/1d/heatcontroldeltax.h \
    parabolic/2d/heatcontrol2d.h \
    parabolic/2d/heatcontrol2delta.h \
    parabolic/2d/heatcontrol2deltaf.h \
    parabolic/2d/heatcontrol2deltax.h \
    parabolic/1d/heatcontrol1.h \
    hyperbolic/hyperbolic2dx.h \
    hyperbolic/hyperbolic1dx.h \
    hyperbolic/hyperboliccontrol1d.h \
    hyperbolic/hyperboliccontrol1d2.h \
    hyperbolic/hyperboliccontrol1d3.h \
    hyperbolic/hyperboliccontrol1dt.h \
    hyperbolic/hyperboliccontrol1d4.h
