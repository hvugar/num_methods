TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
TARGET = main
#QMAKE_CXXFLAGS += -O3
#QMAKE_CXXFLAGS -= -O2
#QMAKE_CXXFLAGS_RELEASE = -O3
#QMAKE_CFLAGS += -O3

INCLUDEPATH += ../minimum
INCLUDEPATH += ../cminimum
LIBS += -L../bin -lminimum -lcminimum

DESTDIR += ../bin

#DEFINES += MINIMUM_LIBRARY
#include(../minimum/minimum.pri)

SOURCES += main.cpp \
    rnfunction/rosenbrock.cpp \
    rnfunction/bealesfunction.cpp \
    rnfunction/boothfunction.cpp \
    control/cfunction.cpp \
    control/cfunction1.cpp \
    control/cfunction2.cpp \
    control/cfunction3.cpp \
    parabolic/1d/heatcontrol.cpp \
    parabolic/1d/heatcontroldeltaf.cpp \
    parabolic/2d/heatcontrol2d.cpp \
    parabolic/2d/heatcontrol2delta.cpp \
    parabolic/2d/heatcontrol2deltaf.cpp \
    parabolic/2d/heatcontrol2deltax.cpp \
    hyperbolic/hyperbolic2dx.cpp \
    hyperbolic/hyperbolic1dx.cpp \
    hyperbolic/hyperboliccontrol1d.cpp \
    hyperbolic/hyperboliccontrol1d4.cpp \
    point/pointcontrol.cpp \
    point/pointcontrol1.cpp \
    point/pointcontrol2.cpp \
    point/pointcontrol11.cpp \
    discrete/discreteheat.cpp \
    discrete/discretehyperbolic.cpp \
    discrete/discretehyperbolic1.cpp \
    hyperbolic/hyperboliccontrolx.cpp \
    parabolic/1d/heatcontroldeltax1.cpp

HEADERS += \
    rnfunction/rosenbrock.h \
    rnfunction/bealesfunction.h \
    rnfunction/boothfunction.h \
    control/cfunction.h \
    control/cfunction1.h \
    control/cfunction2.h \
    control/cfunction3.h \
    parabolic/1d/heatcontrol.h \
    parabolic/1d/heatcontroldeltaf.h \
    parabolic/2d/heatcontrol2d.h \
    parabolic/2d/heatcontrol2delta.h \
    parabolic/2d/heatcontrol2deltaf.h \
    parabolic/2d/heatcontrol2deltax.h \
    hyperbolic/hyperbolic2dx.h \
    hyperbolic/hyperbolic1dx.h \
    hyperbolic/hyperboliccontrol1d.h \
    hyperbolic/hyperboliccontrol1d4.h \
    point/pointcontrol.h \
    point/pointcontrol1.h \
    point/pointcontrol2.h \
    point/pointcontrol11.h \
    discrete/discreteheat.h \
    discrete/discretehyperbolic.h \
    discrete/discretehyperbolic1.h \
    hyperbolic/hyperboliccontrolx.h \
    parabolic/1d/heatcontroldeltax1.h
