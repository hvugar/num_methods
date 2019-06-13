TARGET = problem0H
TEMPLATE = lib
CONFIG  += console
CONFIG  -= app_bundle
CONFIG  -= qt
CONFIG  += shared

OBJECTS_DIR = release/.obj
MOC_DIR     = release/.moc

DEFINES += PROBLEM0H_LIBRARY

DESTDIR = ../bin

INCLUDEPATH += ../minimum
LIBS        += -L../bin -lminimum

HEADERS += problem0h_global.h

#HEADERS += problem0h_common.h
#SOURCES += problem0h_common.cpp

HEADERS += problem0h_solver.h
SOURCES += problem0h_solver.cpp

#HEADERS += testwaveequation.h
#SOURCES += testwaveequation.cpp

#HEADERS += problem2h_exporter.h
#SOURCES += problem2h_exporter.cpp

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
            LIBS        += -L../bin -limaging
        }
    }
    CONFIG(debug, debug|release) { }
}

unix {
    CONFIG(release, debug|release) {}
    CONFIG(debug, debug|release) {}
}

#OTHER_FILES += matlab/*
