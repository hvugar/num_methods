TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += dll

DESTDIR = ../main

DEFINES += MINIMUM_LIBRARY

SOURCES += \
#    minimum.cpp \
#    methods.cpp \
#    method_grad.cpp \
#    method_conj.cpp \
#    method_penalty.cpp \
    function.cpp \
    methods1.cpp \
    minimum1.cpp \
    gradient.cpp

HEADERS += \
    methods.h \
    function.h \
    gradient.h \
    global.h

