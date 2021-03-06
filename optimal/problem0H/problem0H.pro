TARGET   = problem0H
TEMPLATE = lib
CONFIG  += console
CONFIG  -= app_bundle
CONFIG  -= qt
CONFIG  += shared
CONFIG  += c++11

DEFINES += PROBLEM0H_LIBRARY
DEFINES += PROBLEM0H_LIBRARY_ML

DESTDIR = ../bin

#DEFINES += USE_LIB_IMAGING
#DEFINES += USE_LIB_TEXT
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

unix {
    CONFIG(release, debug|release) {}
    CONFIG(debug, debug|release) {}
}

macx-clang {
    CONFIG(release, debug|release) {
        OBJECTS_DIR = release/.obj/macx-clang
        MOC_DIR     = release/.moc/macx-clang

        INCLUDEPATH += ../minimum
        LIBS        += -L../bin -lminimum

        contains(DEFINES, USE_LIB_IMAGING) {
            CONFIG      += qt
            QT          += core gui widgets
            INCLUDEPATH += ../imaging
            LIBS        += -L../bin -limaging
        }

        contains(DEFINES, USE_LIB_XLSX_WRITER) {
            INCLUDEPATH += ../../third-party/macx-clang/x86_64/libxlsxwriter/include
            LIBS        += -L../../third-party/macx-clang/x86_64/libxlsxwriter/lib -lxlsxwriter
        }

        contains(DEFINES, USE_LIB_ZLIB) {
            INCLUDEPATH += ../../third-party/macx-clang/x86_64/zlib/include
            LIBS        += -L../../third-party/macx-clang/x86_64/zlib/lib -lzlibstatic
        }

    }
    CONFIG(debug, debug|release) { }
}

HEADERS += problem0h_global.h

#HEADERS += problem0h_common.h
#SOURCES += problem0h_common.cpp

HEADERS += problem0h_solver.h
SOURCES += problem0h_solver.cpp

HEADERS += problem0h_exporter.h
SOURCES += problem0h_exporter.cpp


#OTHER_FILES += matlab/*
