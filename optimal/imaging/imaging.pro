TARGET   = imaging
TEMPLATE = lib
CONFIG  += shared
CONFIG  -= app_bundle
CONFIG  += c++11

QT     += core gui widgets

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

DESTDIR += ../bin

SOURCES += imaging.cpp

HEADERS += imaging.h
