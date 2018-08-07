TARGET = imaging
TEMPLATE = lib
CONFIG += shared

QT     += core gui widgets

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

DESTDIR += ../bin

SOURCES += imaging.cpp

HEADERS += imaging.h
