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

SOURCES += \
    cmatrix.c \
    function.cpp \
    r1minimize.cpp \
    gradient.cpp \
    gradient_cjt.cpp \
    gradient_sd.cpp \
    penalty.cpp \
#    doublevector.cpp \
    printer.cpp \
    gradient_cs.cpp \
    projection.cpp \
    rungekutta.cpp \
#    tomasmethod.cpp \
    integral.cpp \
    hyperbolicequation.cpp \
    parabolicequation.cpp \
    exceptions.cpp \
    ode1storder.cpp \
#    matrix.cpp \
    matrix2d.cpp \
    vector2d.cpp \
    matrix3d.cpp \
    cmethods.c \
    grid/bvp.cpp \
    grid/ibvp.cpp \
    iibvp.cpp \
    grid/lbvpode.cpp \
    grid/nlbvpode.cpp \
    grid/pibvp.cpp \
    grid/hibvp.cpp \
    grid/hpibvp.cpp \
    grid/grid.cpp \
    grid/bpibvp.cpp

HEADERS += global.h \
    cmethods.h \
    cmatrix.h \
    function.h \
    r1minimize.h \
    gradient.h \
    gradient_cjt.h \
    gradient_sd.h \
    penalty.h \
#    doublevector.h \
    printer.h \
    gradient_cs.h \
    projection.h \
    rungekutta.h \
#    tomasmethod.h \
    integral.h \
    hyperbolicequation.h \
    parabolicequation.h \
    exceptions.h \
    ode1storder.h \
#    matrix.h \
    matrix2d.h \
    vector2d.h \
    matrix3d.h \
    grid/bvp.h \
    grid/ibvp.h \
    iibvp.h \
    grid/lbvpode.h \
    grid/nlbvpode.h \
    grid/pibvp.h \
    grid/hibvp.h \
    grid/hpibvp.h \
    grid/grid.h \
    grid/bpibvp.h
