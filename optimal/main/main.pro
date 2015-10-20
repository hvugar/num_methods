TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH = ../minimum
LIBS += -L../bin -lminimum
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
    heat/1d/heatcontrol.cpp \
    heat/1d/heatcontroldeltaf.cpp \
    heat/1d/heatcontroldeltax.cpp \
    heat/2d/heatcontrol2d.cpp \
    heat/2d/heatcontrol2delta.cpp \
    heat/2d/heatcontrol2deltaf.cpp \
    heat/2d/heatcontrol2deltax.cpp

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
    heat/1d/heatcontrol.h \
    heat/1d/heatcontroldeltaf.h \
    heat/1d/heatcontroldeltax.h \
    heat/2d/heatcontrol2d.h \
    heat/2d/heatcontrol2delta.h \
    heat/2d/heatcontrol2deltaf.h \
    heat/2d/heatcontrol2deltax.h
