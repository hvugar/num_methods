TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += ../minimum
INCLUDEPATH += ../cminimum
LIBS += -L../bin -lminimum -lcminimum
DESTDIR += ../bin

SOURCES += main.cpp \
    example1.cpp \
    example2.cpp \
    example3.cpp \
    example4.cpp \
    example331.cpp

HEADERS += \
    example1.h \
    example2.h \
    example3.h \
    example4.h \
    example331.h
