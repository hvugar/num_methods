TARGET = problem1H
TEMPLATE = lib
CONFIG  += console
CONFIG  -= app_bundle
CONFIG  -= qt
CONFIG  += shared
CONFIG   += c++11

DEFINES += PROBLEM1H_LIBRARY

DESTDIR  = ../bin

#DEFINES += USE_LIB_IMAGING
#DEFINES += USE_LIB_XLSX_WRITER
#DEFINES += USE_LIB_ZLIB

win32-g++ {
    CONFIG(release, debug|release) {
        contains(DEFINES, USE_LIB_XLSX_WRITER) {
            INCLUDEPATH += ../../third-party/MinGW/x86_64/libxlsxwriter/include
            LIBS        += -L../../third-party/MinGW/x86_64/libxlsxwriter/lib -lxlsxwriter
        }
        contains(DEFINES, USE_LIB_ZLIB) {
            INCLUDEPATH += ../../third-party/MinGW/x86_64/zlib/include
            LIBS        += -L../../third-party/MinGW/x86_64/zlib/lib -lzlibstatic
        }
        contains(DEFINES, USE_LIB_IMAGING) {
            CONFIG      += qt
            QT          += core gui widgets
            INCLUDEPATH += ../imaging
            LIBS        += -L../bin -limaging
        }
        OBJECTS_DIR = release/.obj/win32-gcc
        MOC_DIR     = release/.moc/win32-gcc

        INCLUDEPATH += ../minimum
        LIBS        += -L../bin -lminimum

    }
    CONFIG(debug, debug|release) { }
}

win32-msvc* {
    CONFIG(release, debug|release) {
        contains(DEFINES, USE_LIB_XLSX_WRITER) {
            INCLUDEPATH += ../../third-party/VC2015/x86_64/libxlsxwriter/include
            LIBS        += ../../third-party/VC2015/x86_64/libxlsxwriter/lib/xlsxwriter.lib
        }
        contains(DEFINES, USE_LIB_ZLIB) {
            INCLUDEPATH += ../../third-party/VC2015/x86_64/zlib/include
            LIBS        += ../../third-party/VC2015/x86_64/zlib/lib/zlibstatic.lib
        }
        contains(DEFINES, USE_LIB_IMAGING) {
            CONFIG      += qt
            QT          += core gui widgets
            INCLUDEPATH += ../imaging
            LIBS        += ../bin/imaging.lib
        }
        OBJECTS_DIR = release/.obj/win32-msvc
        MOC_DIR     = release/.moc/win32-msvc

        INCLUDEPATH += ../minimum
        LIBS        += ../bin/minimum.lib
    }
    CONFIG(debug, debug|release) { }
}

HEADERS += problem1h_global.h

HEADERS += problem1h_solver.h
SOURCES += problem1h_solver.cpp


#OBJECTS_DIR = release/.obj
#MOC_DIR = release/.moc

#HEADERS += problem1h_global.h \
#           problem1h_common.h \
##           problem2h_solver.h \
#           problem1h_solver1.h \
##           problem2h_solver4.h \
#           problem1h_example.h \
##           problem2h_ibvp.h \
#           problem1h_solver_base.h

#SOURCES += problem1h_common.cpp \
##           problem1h_solver.cpp \
#           problem1h_solver1.cpp \
##           problem2h_solver4.cpp \
#           problem1h_example.cpp \
##           problem2h_ibvp.cpp \
#           problem1h_solver_base.cpp

#HEADERS += problem1h_exporter.h
#SOURCES += problem1h_exporter.cpp

#defined(OLD_SOURCES) {
#    message("Using old sources")
#HEADERS += problem2hnm.h
#SOURCES += problem2hnm.cpp
#}

#OTHER_FILES += matlab/*
