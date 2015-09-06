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
    heatcontrol.cpp \
    pointcontrol.cpp \
    pointcontrol1.cpp \
    pointcontrol2.cpp

HEADERS += \
    utils.h \
    rosenbrock.h \
    bealesfunction.h \
    boothfunction.h \
    cfunction.h \
    cfunction1.h \
    cfunction2.h \
    cfunction3.h \
    heatcontrol.h \
    pointcontrol.h \
    pointcontrol1.h \
    pointcontrol2.h
