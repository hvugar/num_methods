TEMPLATE = lib
TARGET = problem1H
CONFIG  -= app_bundle
#CONFIG  -= qt
CONFIG  += console
CONFIG  += shared
DESTDIR  = ../bin

DEFINES += PROBLEM1H_LIBRARY
DEFINES += USE_IMAGING
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

HEADERS += problem1h_global.h \
           problem1h_common.h \
#           problem2h_solver.h \
           problem1h_solver1.h \
#           problem2h_solver4.h \
           problem1h_example.h \
#           problem2h_ibvp.h \
           problem1h_solver_base.h

SOURCES += problem1h_common.cpp \
#           problem1h_solver.cpp \
           problem1h_solver1.cpp \
#           problem2h_solver4.cpp \
           problem1h_example.cpp \
#           problem2h_ibvp.cpp \
           problem1h_solver_base.cpp

#HEADERS += problem1h_exporter.h
#SOURCES += problem1h_exporter.cpp

defined(OLD_SOURCES) {
    message("Using old sources")
#HEADERS += problem2hnm.h
#SOURCES += problem2hnm.cpp
}

#OTHER_FILES += matlab/*
