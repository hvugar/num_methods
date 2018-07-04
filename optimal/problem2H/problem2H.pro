#-------------------------------------------------
#
# Project created by QtCreator 2018-07-04T14:44:06
#
#-------------------------------------------------

TEMPLATE = lib
TARGET = problem2H
CONFIG  -= app_bundle
CONFIG  -= qt
CONFIG  += console
CONFIG  += dll
DESTDIR  = ../bin

DEFINES += PROBLEM2H_LIBRARY

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

HEADERS += \
        problem2h_global.h \
        problem2h.h

SOURCES += problem2h.cpp
