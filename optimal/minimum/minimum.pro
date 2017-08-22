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
QMAKE_CXXFLAGS += -D__USE_MINGW_ANSI_STDIO

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
    printer.cpp \
    gradient_cs.cpp \
    projection.cpp \
    rungekutta.cpp \
    integral.cpp \
    hyperbolicequation.cpp \
    parabolicequation.cpp \
    exceptions.cpp \
    matrix2d.cpp \
    vector2d.cpp \
    matrix3d.cpp \
    cmethods.c \
    grid/bvp.cpp \
    grid/ibvp.cpp \
    grid/lbvpode.cpp \
    grid/nlbvpode.cpp \
    grid/pibvp.cpp \
    grid/hibvp.cpp \
    grid/hpibvp.cpp \
    grid/grid.cpp \
    grid/bpibvp.cpp \
    grid/nhpibvp.cpp \
    iibvp.cpp \
    gradient/igradient.cpp \
    diffequ.cpp \
    grid/integral1.cpp \
    grid/uniformgrid.cpp \
    nonuniformgrid.cpp \
    ode/cauchyp.cpp \
    utils/matrix.cpp \
    utils/vector.cpp \
    load_sys/islodenlcsm.cpp \
    load_sys/islodenlcsv.cpp \
    load_sys/islodenlcsv2.cpp \
    load_sys/islodenlcs.cpp

HEADERS += global.h \
    cmethods.h \
    cmatrix.h \
    function.h \
    r1minimize.h \
    gradient.h \
    gradient_cjt.h \
    gradient_sd.h \
    penalty.h \
    printer.h \
    gradient_cs.h \
    projection.h \
    rungekutta.h \
    integral.h \
    hyperbolicequation.h \
    parabolicequation.h \
    exceptions.h \
    matrix2d.h \
    vector2d.h \
    matrix3d.h \
    grid/bvp.h \
    grid/ibvp.h \
    grid/lbvpode.h \
    grid/nlbvpode.h \
    grid/pibvp.h \
    grid/hibvp.h \
    grid/hpibvp.h \
    grid/grid.h \
    grid/bpibvp.h \
    grid/nhpibvp.h \
    iibvp.h \
    gradient/igradient.h \
    diffequ.h \
    grid/integral1.h \
    grid/uniformgrid.h \
    nonuniformgrid.h \
    ode/cauchyp.h \
    utils/matrix.h \
    utils/vector.h \
    load_sys/islodenlcsm.h \
    load_sys/islodenlcsv.h \
    load_sys/islodenlcsv2.h \
    load_sys/islodenlcs.h
