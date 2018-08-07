#-------------------------------------------------
#
# Project created by QtCreator 2018-07-04T14:44:06
#
#-------------------------------------------------

TEMPLATE = lib
TARGET = problem2H
CONFIG  -= app_bundle
#CONFIG  -= qt
CONFIG  += console
CONFIG  += shared
DESTDIR  = ../bin
QT       += core gui widgets

DEFINES += PROBLEM2H_LIBRARY

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

INCLUDEPATH += ../imaging
LIBS += -L../bin -limaging

HEADERS += problem2h_global.h common.h

#HEADERS += exporter.h
#SOURCES += exporter.cpp

HEADERS += problem2h.h   problem2hN.h   problem2hnm.h
SOURCES += problem2h.cpp problem2hN.cpp problem2hnm.cpp
