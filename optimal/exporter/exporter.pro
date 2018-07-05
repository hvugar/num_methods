#-------------------------------------------------
#
# Project created by QtCreator 2018-07-04T15:06:29
#
#-------------------------------------------------

TEMPLATE = lib
TARGET = exporter
CONFIG  -= app_bundle
CONFIG  -= qt
CONFIG  += console
CONFIG  += dll
DESTDIR  = ../bin

INCLUDEPATH += ../problem2H
LIBS += -L../bin -lproblem2H

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

DEFINES += EXPORTER_LIBRARY

SOURCES += exporter.cpp

HEADERS += exporter.h
