TARGET = problem0H
TEMPLATE = lib
CONFIG  += console
CONFIG  -= app_bundle
#CONFIG  -= qt
CONFIG  += shared

OBJECTS_DIR = release/.obj
MOC_DIR     = release/.moc

DEFINES += PROBLEM0H_LIBRARY

DESTDIR = ../bin

INCLUDEPATH += ../minimum
LIBS        += -L../bin -lminimum

DEFINES     += USE_IMAGING
INCLUDEPATH += ../imaging
LIBS        += -L../bin -limaging
QT          += core gui widgets

HEADERS += problem0h_global.h

#HEADERS += problem0h_common.h
#SOURCES += problem0h_common.cpp

HEADERS += problem0h_solver.h
SOURCES += problem0h_solver.cpp

#HEADERS += testwaveequation.h
#SOURCES += testwaveequation.cpp

#HEADERS += problem2h_exporter.h
#SOURCES += problem2h_exporter.cpp

DEFINES += USE_LIB_XLSX_WRITER
#defined(USE_LIB_XLSX_WRITER) {
INCLUDEPATH += ../../third-party/VC2015/x86_64/libxlsxwriter/include
#LIBS        += -L../../third-party/VC2015/x86_64/libxlsxwriter/lib -lxlsxwriter
LIBS        += ../../third-party/VC2015/x86_64/libxlsxwriter/lib/xlsxwriter.lib

INCLUDEPATH += ../../third-party/VC2015/x86_64/zlib/include
#LIBS        += -L../../third-party/VC2015/x86_64/zlib/lib -lzlibstatic
LIBS        += ../../third-party/VC2015/x86_64/zlib/lib/zlibstatic.lib
#}

#OTHER_FILES += matlab/*
