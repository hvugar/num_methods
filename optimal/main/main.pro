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
#QMAKE_CXXFLAGS += -Wno-unused-variable
#QMAKE_CXXFLAGS += -O0
#QMAKE_CXXFLAGS += -ffloat-store
#QMAKE_CXXFLAGS += -ffast-math
#QMAKE_CXXFLAGS += -D__USE_MINGW_ANSI_STDIO

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

#INCLUDEPATH += ../imaging
#LIBS += -L../bin -limaging

INCLUDEPATH += ../problem2P
LIBS += -L../bin -lproblem2P

INCLUDEPATH += ../problem2H
LIBS += -L../bin -lproblem2H

INCLUDEPATH += ../problem1H
LIBS += -L../bin -lproblem1H

DESTDIR += ../bin

#include(problem1/problem1.pri)
#include(problem2P/problem2P.pri)
#include(problem4/problem4.pri)
#include(problem5/problem5.pri)
#include(load_sys/load_sys.pri)

SOURCES += main.cpp \
    linearode1storderex1.cpp

#SOURCES += \
#    ivp/nlode1oex1.cpp \
#    nonlinearequationex1.cpp \
#    loadedlinearode1order.cpp \
#    heatequationibvp1.cpp \
#    templatea.cpp \
#    templateb.cpp

#HEADERS += \
#    ivp/nlode1oex1.h \
#    nonlinearequationex1.h \
#    loadedlinearode1order.h \
#    heatequationibvp1.h \
#    templatea.h \
#    templateb.h

OBJECTS_DIR = release/.obj
MOC_DIR = release/.moc

HEADERS += \
    linearode1storderex1.h
