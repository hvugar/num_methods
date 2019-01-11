TEMPLATE = lib
TARGET = problem2H
CONFIG  -= app_bundle
#CONFIG  -= qt
CONFIG  += console
CONFIG  += shared
DESTDIR  = ../bin

DEFINES += PROBLEM2H_LIBRARY
#DEFINES += USE_IMAGING
#DEFINES += OLD_SOURCES

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

#defined(USE_IMAGING) {
#message("Using imaging library")
INCLUDEPATH += ../imaging
LIBS        += -L../bin -limaging
QT          += core gui widgets
#}

OBJECTS_DIR = release/.obj
MOC_DIR = release/.moc

HEADERS += problem2h_global.h \
           problem2h_common.h \
           problem2h_solver.h \
           problem2h_solver1.h \
           problem2h_solver4.h \
           problem2h_example.h \
           problem2h_ibvp.h \
    problem2h_solver_base.h

SOURCES += problem2h_common.cpp \
           problem2h_solver.cpp \
           problem2h_solver1.cpp \
           problem2h_solver4.cpp \
           problem2h_example.cpp \
           problem2h_ibvp.cpp \
    problem2h_solver_base.cpp

HEADERS += problem2h_exporter.h
SOURCES += problem2h_exporter.cpp

defined(OLD_SOURCES) {
    message("Using old sources")
#HEADERS += problem2hnm.h
#SOURCES += problem2hnm.cpp
}



OTHER_FILES += matlab/*
