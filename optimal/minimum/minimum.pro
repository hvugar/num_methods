TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += dll
#DEFINES += OS_UNIX

#QMAKE_CXXFLAGS += -O3
#QMAKE_CXXFLAGS -= -O2
#QMAKE_CXXFLAGS_RELEASE = -O3
#QMAKE_CFLAGS += -O3

DEFINES += MINIMUM_LIBRARY

DESTDIR += ../bin

INCLUDEPATH += ../cminimum
LIBS += -L../bin -lcminimum

SOURCES += \
    function.cpp \
    r1minimize.cpp \
    gradient.cpp \
    gradient_cjt.cpp \
    gradient_sd.cpp \
    penalty.cpp \
    doublevector.cpp \
    printer.cpp \
    gradient_cs.cpp \
    projection.cpp \
    rungekutta.cpp \
    tomasmethod.cpp \
    integral.cpp \
    hyperbolicequation.cpp \
    parabolicequation.cpp \
    exceptions.cpp

HEADERS += global.h \
    function.h \
    r1minimize.h \
    gradient.h \
    gradient_cjt.h \
    gradient_sd.h \
    penalty.h \
    doublevector.h \
    printer.h \
    gradient_cs.h \
    projection.h \
    rungekutta.h \
    tomasmethod.h \
    integral.h \
    hyperbolicequation.h \
    parabolicequation.h \
    exceptions.h
