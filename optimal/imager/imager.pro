#-------------------------------------------------
#
# Project created by QtCreator 2015-10-20T09:47:22
#
#-------------------------------------------------

QT       += core
QT       += gui

TARGET = imager
CONFIG   += console
CONFIG   -= app_bundle
DESTDIR += ../bin
TEMPLATE = app

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum -lcminimum


SOURCES += main.cpp widget2.cpp

HEADERS += \
    widget2.h
