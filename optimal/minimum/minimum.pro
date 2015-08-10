TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += dll

DEFINES += MINIMUM_LIBRARY

DESTDIR += ../bin

SOURCES += \
    function.cpp \
    r1minimize.cpp \
    gradient.cpp \
    gradient_cjt.cpp \
    gradient_sd.cpp \
    gradient_prj.cpp \
    penalty.cpp \
    doublevector.cpp \
    gridmethod.cpp \
    printer.cpp \
    gradient_cs.cpp

HEADERS += global.h \
    function.h \
    r1minimize.h \
    gradient.h \
    gradient_cjt.h \
    gradient_sd.h \
    gradient_prj.h \
    penalty.h \
    doublevector.h \
    gridmethod.h \
    printer.h \
    gradient_cs.h

DISTFILES += minimum.pri
