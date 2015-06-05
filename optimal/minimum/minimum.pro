TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += dll

DEFINES += MINIMUM_LIBRARY

DESTDIR += ../bin
DESTDIR_TARGET += ../lib

SOURCES += \
    function.cpp \
    methods1.cpp \
    minimum1.cpp \
    r1minimize.cpp \
    gradient.cpp \
    fpgradient.cpp

HEADERS += global.h \
    methods.h \
    function.h \
    r1minimize.h \
    gradient.h \
    fpgradient.h
