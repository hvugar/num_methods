TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH = ../minimum
LIBS += -L. -l../bin/minimum

DESTDIR += ../bin

SOURCES += main.cpp \
    samplecontrol.cpp \
    cfunction1.cpp \
    cfunction2.cpp \
    utils.cpp

HEADERS += \
    samplecontrol.h \
    cfunction1.h \
    cfunction2.h \
    utils.h

