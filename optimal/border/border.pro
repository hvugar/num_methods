TARGET   = border
TEMPLATE = lib
CONFIG  += console
CONFIG  -= app_bundle
CONFIG  -= qt
CONFIG  += shared

DEFINES += BORDER_LIBRARY

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

unix {
    CONFIG(release, debug|release) {}
    CONFIG(debug, debug|release) {}
}

HEADERS += border_global.h

HEADERS += \
           waveequationibvp.h \
           newtonheatequationex1.h

SOURCES += \
           waveequationibvp.cpp \
           newtonheatequationex1.cpp

#HEADERS += \
#           ode/firstordernonlinearodeex1.h \
#           ode/firstorderlinearodeex1.h \
#           ode/secondorderlinearodeex1.h

#SOURCES += \
#           ode/firstordernonlinearodeex1.cpp \
#           ode/firstorderlinearodeex1.cpp \
#           ode/secondorderlinearodeex1.cpp

# Input
HEADERS += \
#           borderparabolicd.h \
#           borderparabolicn.h \
#           borderparabolic2d.h \
#           borderhyperbolic2d.h \
#           borderhyperbolic2d1.h \
#           grid/parabolicibvp1.h \

SOURCES +=
#           borderparabolicd.cpp \
#           borderparabolicn.cpp \
#           borderparabolic2d.cpp \
#           borderhyperbolic2d.cpp \
#           borderhyperbolic2d1.cpp \
#           grid/parabolicibvp1.cpp \
