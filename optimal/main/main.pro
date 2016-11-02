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
#INCLUDEPATH += ../imaging
#LIBS += -L../bin -limaging

DESTDIR += ../bin

QT += core gui widgets

SOURCES += main.cpp \
    loadedsystems.cpp \
    example1.cpp \
    example2.cpp \
    problem1.cpp \
    bordertest.cpp

HEADERS += \
    loadedsystems.h \
    example1.h \
    example2.h \
    problem1.h \
    bordertest.h
