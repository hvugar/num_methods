TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
#CONFIG -= qt
TARGET = main
#QMAKE_CXXFLAGS += -O3
#QMAKE_CXXFLAGS -= -O2
#QMAKE_CXXFLAGS_RELEASE = -O3
#QMAKE_CFLAGS += -O3
#QMAKE_CXXFLAGS += -Wunused
#QMAKE_CXXFLAGS += -Werror
#QMAKE_CXXFLAGS -= -Wall
QMAKE_CXXFLAGS += -Wno-unused-variable

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum
INCLUDEPATH += ../imaging
LIBS += -L../bin -limaging

DESTDIR += ../bin

QT += core gui widgets

SOURCES += main.cpp \
    problem1.cpp \
    problem3.cpp \
    problem1k.cpp \
    problem1z.cpp \
    problem1x.cpp \
    problem1x1.cpp \
    problem1x2.cpp \
    problem1kz.cpp \
    loadedsystems.cpp \
    r1functionvector.cpp \
    example1.cpp

HEADERS += \
    problem1.h \
    problem3.h \
    problem1k.h \
    problem1x.h \
    problem1x1.h \
    problem1x2.h \
    problem1kz.h \
    problem1z.h \
    loadedsystems.h \
    r1functionvector.h \
    example1.h
