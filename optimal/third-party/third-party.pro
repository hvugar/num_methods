TEMPLATE = lib
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
CONFIG  += dll
DESTDIR  = ../bin

SOURCES +=

INCLUDEPATH += ../../../../../third-party/jmcnamara/libxlsxwriter/include
LIBS +=      -L../../../../../third-party/jmcnamara/libxlsxwriter/lib -lxlsxwriter
LIBS += -L. -lz
