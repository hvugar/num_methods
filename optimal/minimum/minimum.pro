TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += dll

DEFINES += MINIMUM_LIBRARY

DESTDIR += ../bin

SOURCES += \
    function.cpp \
    methods1.cpp \
    minimum1.cpp \
    r1minimize.cpp \
    gradient.cpp \
    cjtgradient.cpp \
    prjgradient.cpp \
    sdgradient.cpp

HEADERS += global.h \
    methods.h \
    function.h \
    r1minimize.h \
    gradient.h \
    cjtgradient.h \
    prjgradient.h \
    sdgradient.h
