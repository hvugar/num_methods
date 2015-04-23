TEMPLATE = lib
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += dll

DESTDIR = ../main

SOURCES += minimum.cpp \
    methods.cpp

HEADERS += \
    methods.h

