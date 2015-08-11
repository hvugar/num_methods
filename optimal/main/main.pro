TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH = ../minimum
LIBS += -L../bin -lminimum
DESTDIR += ../bin

#include(../minimum/minimum.pri)

SOURCES += main.cpp \
    cfunction1.cpp \
    cfunction2.cpp \
    utils.cpp \
#    heatcontrol.cpp \
    rosenbrock.cpp \
    cfunction.cpp \
    bealesfunction.cpp

HEADERS += \
    cfunction1.h \
    cfunction2.h \
    utils.h \
#    heatcontrol.h \
    rosenbrock.h \
    cfunction.h \
    bealesfunction.h

