TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += ../minimum
INCLUDEPATH += ../cminimum
LIBS += -L../bin -lminimum -lcminimum
DESTDIR += ../bin

SOURCES += main.cpp \
    borderparabolic1d.cpp

HEADERS += \
    borderparabolic1d.h
