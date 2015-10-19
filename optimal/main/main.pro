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
    heat/heatcontrol.cpp \
    heat/heatcontrol2d.cpp \
    heat/heatcontrol2delta.cpp \
    heat/heatcontrol2deltaf.cpp \
    heat/heatcontroldeltaf.cpp \
    heat/heatcontroldeltax.cpp \
    heat/heatcontrol2deltax.cpp

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
    heat/heatcontrol.h \
    heat/heatcontrol2d.h \
    heat/heatcontrol2delta.h \
    heat/heatcontrol2deltaf.h \
    heat/heatcontroldeltaf.h \
    heat/heatcontroldeltax.h \
    heat/heatcontrol2deltax.h
