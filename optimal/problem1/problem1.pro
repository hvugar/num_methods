TEMPLATE = lib
TARGET   = problem1
CONFIG  -= app_bundle
#CONFIG -= qt
CONFIG  += console
CONFIG  += shared
DESTDIR  = ../bin
QT      += core gui widgets

DEFINES += PROBLEM1H_LIBRARY

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

INCLUDEPATH += ../imaging
LIBS += -L../bin -limaging

SOURCES += \
    art_problem1.cpp \
    ibackwardloadedheateqauation.cpp \
    iloadedheatequation.cpp \
    iproblem1.cpp \
    loadedheatequation.cpp \
    problem1L1.cpp \
    problem1L2.cpp \
    problem1L3.cpp

HEADERS += \
    art_problem1.h \
    ibackwardloadedheateqauation.h \
    iloadedheatequation.h \
    iproblem1.h \
    loadedheatequation.h \
    problem1L1.h \
    problem1L2.h \
    problem1L3.h
