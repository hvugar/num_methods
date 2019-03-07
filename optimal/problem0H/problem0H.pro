TEMPLATE = lib
TARGET = problem0H
CONFIG  -= app_bundle
#CONFIG  -= qt
CONFIG  += console
CONFIG  += shared
DESTDIR  = ../bin

DEFINES += PROBLEM2H_LIBRARY
DEFINES += USE_IMAGING

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

HEADERS += problem0h_global.h

#HEADERS += problem0h_common.h
#SOURCES += problem0h_common.cpp

HEADERS += problem0h_solver.h
SOURCES += problem0h_solver.cpp

#HEADERS += problem2h_exporter.h
#SOURCES += problem2h_exporter.cpp

DEFINES += USE_LIB_XLSX_WRITER
#defined(USE_LIB_XLSX_WRITER) {
INCLUDEPATH += ../../third-party/MinGW/x86_64/libxlsxwriter/include
LIBS        += -L../../third-party/MinGW/x86_64/libxlsxwriter/lib -lxlsxwriter

INCLUDEPATH += ../../third-party/MinGW/x86_64/zlib/include
LIBS        += -L../../third-party/MinGW/x86_64/zlib/lib -lzlibstatic
#}

#OTHER_FILES += matlab/*
