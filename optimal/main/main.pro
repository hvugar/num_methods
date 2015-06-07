TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH = ../minimum
LIBS += -L. -l../bin/minimum

DESTDIR += ../bin

SOURCES += main.cpp \
    samplecontrol.cpp \
    cfunction.cpp

HEADERS += \
    samplecontrol.h \
    cfunction.h

