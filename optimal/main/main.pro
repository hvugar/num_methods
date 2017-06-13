TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
TARGET = main
#QMAKE_CXXFLAGS += -O2
#QMAKE_CXXFLAGS_RELEASE = -O3
#QMAKE_CFLAGS += -O3
#QMAKE_CXXFLAGS += -Wunused
#QMAKE_CXXFLAGS += -Werror
#QMAKE_CXXFLAGS -= -Wall
QMAKE_CXXFLAGS += -Wno-unused-variable
#QMAKE_CXXFLAGS += -O0
#QMAKE_CXXFLAGS += -ffloat-store
#QMAKE_CXXFLAGS += -ffast-math

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

INCLUDEPATH += ../border
LIBS += -L../bin -lborder

#INCLUDEPATH += ../hyperbolic
#LIBS += -L../bin -lhyperbolic

INCLUDEPATH += ../parabolic
LIBS += -L../bin -lparabolic

INCLUDEPATH += ../rnfunction
LIBS += -L../bin -lrnfunction

#INCLUDEPATH += ../imaging
#LIBS += -L../bin -limaging

DESTDIR += ../bin

SOURCES += main.cpp \
    loadedsystems.cpp \
    example1.cpp \
    example2.cpp \
    bordertest.cpp \
    bordertest1.cpp \
    example3.cpp \
    example4.cpp \
    example5.cpp \
    sampleboundaryproblem1.cpp \
    example6.cpp \
    example7.cpp \
    problem1/iproblem1.cpp \
    problem1/problem1L2.cpp \
    problem1/problem1L3.cpp \
    problem1/art_problem1.cpp \
    problem1/problem1L1.cpp \
    high_order/singledifequ.cpp \
    high_order/systemdifequ.cpp

HEADERS += \
    loadedsystems.h \
    example1.h \
    example2.h \
    bordertest.h \
    bordertest1.h \
    example3.h \
    example4.h \
    example5.h \
    sampleboundaryproblem1.h \
    example6.h \
    example7.h \
    problem1/iproblem1.h \
    problem1/problem1L2.h \
    problem1/problem1L3.h \
    problem1/art_problem1.h \
    problem1/problem1L1.h \
    high_order/singledifequ.h \
    high_order/systemdifequ.h
