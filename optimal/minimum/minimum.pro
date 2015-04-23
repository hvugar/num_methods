TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += dll

DESTDIR = ../main

SOURCES += minimum.cpp \
    methods.cpp \
    method_grad.cpp \
    method_conj.cpp \
    method_penalty.cpp

HEADERS += \
    methods.h

