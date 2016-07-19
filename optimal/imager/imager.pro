#-------------------------------------------------
#
# Project created by QtCreator 2015-10-20T09:47:22
#
#-------------------------------------------------

QT       += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = imager
CONFIG   += console
CONFIG   -= app_bundle
DESTDIR += ../bin
TEMPLATE = app

INCLUDEPATH += ../minimum
INCLUDEPATH += ../cminimum
LIBS += -L../bin -lminimum -lcminimum

SOURCES += main.cpp widget2.cpp sample1.cpp heatimager.cpp

HEADERS += widget2.h sample1.h heatimager.h
