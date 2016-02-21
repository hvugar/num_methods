TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
TARGET = main
#QMAKE_CXXFLAGS += -O3
#QMAKE_CXXFLAGS -= -O2
#QMAKE_CXXFLAGS_RELEASE = -O3
#QMAKE_CFLAGS += -O3
#QMAKE_CXXFLAGS += -Wunused
#QMAKE_CXXFLAGS += -Werror
#QMAKE_CXXFLAGS -= -Wall

INCLUDEPATH += ../minimum
INCLUDEPATH += ../cminimum
LIBS += -L../bin -lminimum -lcminimum

INCLUDEPATH += D:/libxl-3.6.5.0/include_cpp
LIBS += D:/libxl-3.6.5.0/lib/libxl.lib

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
    parabolic/1d/heatcontroldeltax.cpp \
    parabolic/2d/heatcontrol2d.cpp \
    parabolic/2d/heatcontrol2delta.cpp \
    parabolic/2d/heatcontrol2deltaf.cpp \
    parabolic/2d/heatcontrol2deltax.cpp \
    hyperbolic/1d/hyperbolic1dx.cpp \
    hyperbolic/1d/hyperboliccontrol1d.cpp \
    hyperbolic/1d/hyperboliccontrolh.cpp \
    hyperbolic/1d/hyperboliccontrolx.cpp \
    hyperbolic/2d/hyperboliccontrol2d.cpp \
    hyperbolic/2d/hyperboliccontrol2dm.cpp \
    hyperbolic/2d/hyperboliccontrol2dmx.cpp \
    hyperbolic/2d/hyperboliccontrol2dmv.cpp \
    point/pointcontrol.cpp \
    point/pointcontrol1.cpp \
    point/pointcontrol2.cpp \
    point/pointcontrol11.cpp \
    discrete/discreteheat.cpp \
    discrete/discretehyperbolic.cpp \
    discrete/discretehyperbolic1.cpp \
    border/borderparabolic2d.cpp \
    border/borderhyperbolic2d.cpp

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
    parabolic/1d/heatcontroldeltax.h \
    parabolic/2d/heatcontrol2d.h \
    parabolic/2d/heatcontrol2delta.h \
    parabolic/2d/heatcontrol2deltaf.h \
    parabolic/2d/heatcontrol2deltax.h \
    hyperbolic/1d/hyperbolic1dx.h \
    hyperbolic/1d/hyperboliccontrol1d.h \
    hyperbolic/1d/hyperboliccontrolh.h \
    hyperbolic/1d/hyperboliccontrolx.h \
    hyperbolic/2d/hyperboliccontrol2d.h \
    hyperbolic/2d/hyperboliccontrol2dm.h \
    hyperbolic/2d/hyperboliccontrol2dmx.h \
    hyperbolic/2d/hyperboliccontrol2dmv.h \
    point/pointcontrol.h \
    point/pointcontrol1.h \
    point/pointcontrol2.h \
    point/pointcontrol11.h \
    discrete/discreteheat.h \
    discrete/discretehyperbolic.h \
    discrete/discretehyperbolic1.h \
    border/borderparabolic2d.h \
    border/borderhyperbolic2d.h \
    headers.h
