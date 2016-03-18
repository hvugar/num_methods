TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += ../minimum
INCLUDEPATH += ../cminimum
LIBS += -L../bin -lminimum -lcminimum
DESTDIR += ../bin

SOURCES += main.cpp \
    example311.cpp \
    example312.cpp \
    example321.cpp \
    example322.cpp \
    example331.cpp \
    example332.cpp \
    example333.cpp \
    example334.cpp

HEADERS += \
    example311.h \
    example312.h \
    example321.h \
    example322.h \
    example331.h \
    example332.h \
    example333.h \
    example334.h
