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
    heat/heatcontrol.cpp \
    heat/heatcontroldelta.cpp \
    pointcontrol.cpp \
    pointcontrol1.cpp \
    pointcontrol2.cpp \
    heat/headcontrol2d.cpp \
    heat/headcontrol2d1.cpp \
    heat/heatcontrol2d.cpp

HEADERS += \
    utils.h \
    rosenbrock.h \
    bealesfunction.h \
    boothfunction.h \
    cfunction.h \
    cfunction1.h \
    cfunction2.h \
    cfunction3.h \
    heat/heatcontrol.h \
    heat/heatcontroldelta.h \
    pointcontrol.h \
    pointcontrol1.h \
    pointcontrol2.h \
    heat/headcontrol2d.h \
    heat/headcontrol2d1.h \
    heat/heatcontrol2d.h
