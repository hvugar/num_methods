TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
#CONFIG -= qt
QT       += core gui widgets
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
QMAKE_CXXFLAGS += -D__USE_MINGW_ANSI_STDIO

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

INCLUDEPATH += ../border
LIBS += -L../bin -lborder

#INCLUDEPATH += ../hyperbolic
#LIBS += -L../bin -lhyperbolic

#INCLUDEPATH += ../parabolic
#LIBS += -L../bin -lparabolic

#INCLUDEPATH += ../rnfunction
#LIBS += -L../bin -lrnfunction

INCLUDEPATH += ../imaging
LIBS += -L../bin -limaging

DESTDIR += ../bin

SOURCES += main.cpp \
    load_sys/slodenlcsm.cpp \
    load_sys/slodenlcsv.cpp \
    problem1/iproblem1.cpp \
    problem1/problem1L2.cpp \
    problem1/problem1L3.cpp \
    problem1/art_problem1.cpp \
    problem1/problem1L1.cpp \
    problem1/iloadedheatequation.cpp \
    problem1/loadedheatequation.cpp \
    problem1/ibackwardloadedheateqauation.cpp \
    problem4/problem4ex1.cpp \
    problem4/zetta0.cpp \
    numintegralexp1.cpp \
    load_sys/slodenlcsv2.cpp \
    problem4/zettai.cpp \
    ivp/nlode1oex1.cpp \
    load_sys/lode1oex1.cpp \
    matrixtest.cpp \
    nonlinearequationex1.cpp \
    problem4/problem4ex2.cpp \
    loadedlinearode1order.cpp \
    problem5/nllparabolic.cpp \
    problem2/1d/problem2.cpp \
    problem2/1d/iproblem2forward.cpp \
    problem2/1d/iproblem2backward.cpp \
    problem2/2d/iproblem2forward2d.cpp \
    problem2/2d/iproblem2backward2d.cpp \
    problem2/2d/problem22d.cpp \
    problem2/2d/problem2setting.cpp \
    problem2/2d/abstractproblem22d.cpp \
    problem2/2d/problem22dex1.cpp

HEADERS += \
    load_sys/slodenlcsm.h \
    load_sys/slodenlcsv.h \
    problem1/iproblem1.h \
    problem1/problem1L2.h \
    problem1/problem1L3.h \
    problem1/art_problem1.h \
    problem1/problem1L1.h \
    problem1/iloadedheatequation.h \
    problem1/loadedheatequation.h \
    problem1/ibackwardloadedheateqauation.h \
    problem4/problem4ex1.h \
    problem4/zetta0.h \
    numintegralexp1.h \
    load_sys/slodenlcsv2.h \
    problem4/zettai.h \
    ivp/nlode1oex1.h \
    load_sys/lode1oex1.h \
    matrixtest.h \
    nonlinearequationex1.h \
    problem4/problem4ex2.h \
    loadedlinearode1order.h \
    problem5/nllparabolic.h \
    problem2/1d/problem2.h \
    problem2/1d/iproblem2backward.h \
    problem2/1d/iproblem2forward.h \
    problem2/2d/iproblem2forward2d.h \
    problem2/2d/iproblem2backward2d.h \
    problem2/2d/problem22d.h \
    problem2/2d/problem2setting.h \
    problem2/2d/abstractproblem22d.h \
    problem2/2d/problem22dex1.h
