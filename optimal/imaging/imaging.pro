TARGET = imaging
TEMPLATE = lib
CONFIG += static

QT       += core gui widgets
CONFIG += static

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

DESTDIR += ../bin

SOURCES += imaging.cpp

HEADERS += imaging.h
