#-------------------------------------------------
#
# Project created by QtCreator 2016-08-12T10:39:11
#
#-------------------------------------------------

QT       += core gui widgets

TARGET = imaging
TEMPLATE = lib

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

#include(../minimum/minimum.pri)

DESTDIR += ../bin

SOURCES += imaging.cpp

HEADERS += imaging.h
