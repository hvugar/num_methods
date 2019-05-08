TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

OBJECTS_DIR = release/.obj
MOC_DIR = release/.moc

SOURCES += \
        main.cpp
