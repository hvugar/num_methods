TEMPLATE = app
CONFIG  += console
CONFIG  -= app_bundle
#CONFIG  -= qt
QT      += core gui widgets
TARGET  = main
CONFIG  += c++11

#QMAKE_CXXFLAGS += -O2
#QMAKE_CXXFLAGS_RELEASE = -O3
#QMAKE_CFLAGS += -O3
#QMAKE_CXXFLAGS += -Wunused
#QMAKE_CXXFLAGS += -Werror
#QMAKE_CXXFLAGS -= -Wall
#QMAKE_CXXFLAGS += -Wno-unused-variable
#QMAKE_CXXFLAGS += -O0
#QMAKE_CXXFLAGS += -ffloat-store
#QMAKE_CXXFLAGS += -ffast-math
#QMAKE_CXXFLAGS += -D__USE_MINGW_ANSI_STDIO

DESTDIR = ../bin

INCLUDEPATH += ../border
LIBS += -L../bin -lborder

#INCLUDEPATH += ../hyperbolic
#LIBS += -L../bin -lhyperbolic

#INCLUDEPATH += ../parabolic
#LIBS += -L../bin -lparabolic

#INCLUDEPATH += ../rnfunction
#LIBS += -L../bin -lrnfunction

#INCLUDEPATH += ../imaging
#LIBS += -L../bin -limaging

#INCLUDEPATH += ../problem2P
#LIBS += -L../bin -lproblem2P

#INCLUDEPATH += ../problem1H
#LIBS += -L../bin -lproblem1H

#include(problem4/problem4.pri)
#include(problem5/problem5.pri)
#include(load_sys/load_sys.pri)

SOURCES += main.cpp \
    conjugate_gradinet_test.cpp \
    test/delta_grid_2d_ext1.cpp \
    test/nonlinear_equation_ex1.cpp

SOURCES +=
HEADERS += \
    conjugate_gradinet_test.h \
    test/delta_grid_2d_ext1.h \
    test/nonlinear_equation_ex1.h

SOURCES +=
HEADERS +=

win32-g++ {
    CONFIG(release, debug|release) {
        OBJECTS_DIR = release/.obj/win32-gcc
        MOC_DIR     = release/.moc/win32-gcc

        INCLUDEPATH += ../minimum
        LIBS += -L../bin -lminimum

        INCLUDEPATH += ../problem0H
        LIBS += -L../bin -lproblem0H

        INCLUDEPATH += ../problem2H
        LIBS += -L../bin -lproblem2H
    }
    CONFIG(debug, debug|release) {
        OBJECTS_DIR = debug/.obj/win32-gcc
        MOC_DIR     = debug/.moc/win32-gcc
    }
}

win32-msvc* {
    CONFIG(release, debug|release) {
        OBJECTS_DIR = release/.obj/win32-msvc
        MOC_DIR     = release/.moc/win32-msvc

        INCLUDEPATH += ../minimum
        LIBS += ../bin/minimum.lib

        INCLUDEPATH += ../problem0H
        LIBS += ../bin/problem0H.lib

        INCLUDEPATH += ../problem2H
        LIBS += ../bin/problem2H.lib
    }
    CONFIG(debug, debug|release) {
        OBJECTS_DIR = debug/.obj/win32-msvc
        MOC_DIR     = debug/.moc/win32-msvc
    }
}

unix {
    CONFIG(release, debug|release) {}
    CONFIG(debug, debug|release) {}
}

macx-clang{
    CONFIG(release, debug|release) {
        OBJECTS_DIR = release/.obj/macx-clang
        MOC_DIR     = release/.moc/macx-clang

        INCLUDEPATH += ../minimum
        LIBS += -L../bin -lminimum

#        INCLUDEPATH += ../problem0H
#        LIBS += -L../bin -lproblem0H

#        INCLUDEPATH += ../problem2H
#        LIBS += -L../bin -lproblem2H
    }

    CONFIG(debug, debug|release) {
        OBJECTS_DIR = debug/.obj/macx-clang
        MOC_DIR     = debug/.moc/macx-clang
    }
}
