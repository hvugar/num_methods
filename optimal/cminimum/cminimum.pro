TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += dll

DEFINES += C_MINIMUM_LIBRARY
DESTDIR += ../bin

SOURCES += \
    cmethods.c

HEADERS += \
    cmethods.h \
    cglobal.h

