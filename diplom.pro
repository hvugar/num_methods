TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

DESTDIR += bin

SOURCES += main.cpp \
    example31.cpp \
    example32.cpp \
    example33.cpp \
    example34.cpp \
    example35.cpp \
    example36.cpp \
    example37.cpp \
    example371.cpp

HEADERS += \
    example31.h \
    example32.h \
    example33.h \
    example34.h \
    example35.h \
    example36.h \
    example37.h \
    example371.h

SOURCES += \
    doublevector.cpp \
    parabolicequation.cpp \
    printer.cpp \
    gradient.cpp \
    gradient_cjt.cpp \
    function.cpp \
    r1minimize.cpp

HEADERS += \
    doublevector.h \
    parabolicequation.h \
    printer.h \
    gradient.h \
    gradient_cjt.h \
    function.h \
    r1minimize.h
