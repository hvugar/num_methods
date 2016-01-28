TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += ../minimum
INCLUDEPATH += ../cminimum
LIBS += -L../bin -lminimum -lcminimum

DESTDIR += ../bin

SOURCES += main.cpp \
    heatdeltacenter.cpp

HEADERS += \
    heatdeltacenter.h

