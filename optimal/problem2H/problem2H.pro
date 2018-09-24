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

HEADERS += problem2h_global.h
HEADERS += common.h

HEADERS += problem2h_exporter.h
SOURCES += problem2h_exporter.cpp

HEADERS += problem2h_solver.h problem2h_example.h
SOURCES += problem2h_solver.cpp problem2h_example.cpp

#HEADERS += problem2hnm.h
#SOURCES += problem2hnm.cpp

OTHER_FILES += matlab/*

OBJECTS_DIR = release/.obj
MOC_DIR = release/.moc
