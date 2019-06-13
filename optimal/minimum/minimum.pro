TARGET = minimum
TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += shared

OBJECTS_DIR = release/.obj
MOC_DIR = release/.moc

DEFINES += MINIMUM_LIBRARY

DESTDIR = ../bin

SOURCES += \
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
    #ode/cauchyp.cpp \
    ode/diffequ.cpp \
    pde_old/iibvp.cpp \
    pde_old/hyperbolicequation.cpp \
    pde_old/parabolicequation.cpp \
    exceptions.cpp \
    vector2d.cpp \
    matrix2d.cpp \
    matrix3d.cpp \
    cmethods.c \
    grid/bvp.cpp \
    grid/ibvp.cpp \
    grid/pibvp.cpp \
    grid/hibvp.cpp \
    grid/grid.cpp \
    grid/integral1.cpp \
    nonuniformgrid.cpp \
    utils/matrix.cpp \
    utils/vector.cpp \
    load_sys/islodenlcsm.cpp \
    load_sys/islodenlcsv.cpp \
    load_sys/islodenlcsv2.cpp \
    load_sys/islodenlcs.cpp \
    ode/nlode1o.cpp \
    ode/lode1o.cpp \
    linearequation.cpp \
    nonlinearequation.cpp \
    utils/random.cpp \
    ode/lode2o.cpp \
    ode/nlode2o.cpp \
    vectornormalizer.cpp \
    benchmark.cpp \
    deltagrid.cpp \
    interpolation.cpp \
    grid/ivp.cpp

HEADERS += global.h \
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
    pde_old/iibvp.h \
    pde_old/hyperbolicequation.h \
    pde_old/parabolicequation.h \
    exceptions.h \
    matrix2d.h \
    vector2d.h \
    matrix3d.h \
    cmethods.h \
    grid/bvp.h \
    grid/ibvp.h \
    grid/pibvp.h \
    grid/hibvp.h \
    grid/grid.h \
    grid/integral1.h \
    nonuniformgrid.h \
    #ode/cauchyp.h \
    ode/diffequ.h \
    utils/matrix.h \
    utils/vector.h \
    load_sys/islodenlcsm.h \
    load_sys/islodenlcsv.h \
    load_sys/islodenlcsv2.h \
    load_sys/islodenlcs.h \
    ode/nlode1o.h \
    ode/lode1o.h \
    linearequation.h \
    nonlinearequation.h \
    utils/random.h \
    ode/lode2o.h \
    ode/nlode2o.h \
    vectornormalizer.h \
    benchmark.h \
    deltagrid.h \
    interpolation.h \
    grid/ivp.h
