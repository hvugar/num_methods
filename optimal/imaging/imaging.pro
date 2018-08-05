TARGET = imaging
TEMPLATE = lib
CONFIG += static

QT       += core gui widgets

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

DESTDIR += ../bin

SOURCES += imaging.cpp

HEADERS += imaging.h
