QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = grapher
TEMPLATE = app

DESTDIR += ../bin


SOURCES += main.cpp\
        widget.cpp

HEADERS  += widget.h
