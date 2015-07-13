TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH = ../minimum
DEPENDPATH += ../minimum
LIBS += -L../bin -lminimum
DESTDIR += ../bin

SOURCES += main.cpp \
    samplecontrol.cpp \
    cfunction1.cpp \
    cfunction2.cpp \
    utils.cpp \
    headcontrol.cpp

HEADERS += \
    samplecontrol.h \
    cfunction1.h \
    cfunction2.h \
    utils.h \
    headcontrol.h

