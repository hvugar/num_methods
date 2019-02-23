######################################################################
# Automatically generated by qmake (3.0) Sat Sep 3 16:31:09 2016
######################################################################

TEMPLATE = lib
TARGET   = border
CONFIG  += console
CONFIG  -= app_bundle
#CONFIG  -= qt
CONFIG  += shared
DESTDIR  = ../bin

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

INCLUDEPATH += ../imaging
LIBS        += -L../bin -limaging
QT          += core gui widgets

DEFINES += MINIMUM_LIBRARY

# Input
HEADERS +=
#           borderparabolicd.h \
#           borderparabolicn.h \
#           borderparabolic2d.h \
#           borderhyperbolic2d.h \
#           borderhyperbolic2d1.h \
#           grid/parabolicibvp1.h \

SOURCES +=
#           borderparabolicd.cpp \
#           borderparabolicn.cpp \
#           borderparabolic2d.cpp \
#           borderhyperbolic2d.cpp \
#           borderhyperbolic2d1.cpp \
#           grid/parabolicibvp1.cpp \

HEADERS += \
           grid/hyperbolicibvp1.h \
           grid/newtonheatequationex1.h \
           grid/secondorderlinearodeex1.h

SOURCES += \
           grid/hyperbolicibvp1.cpp \
           grid/newtonheatequationex1.cpp \
           grid/secondorderlinearodeex1.cpp


HEADERS += \
           ode/firstordernonlinearodeex1.h \
           ode/linearode1storderex1.h \
           ode/nonlinearode1storderex1.h

SOURCES += \
           ode/firstordernonlinearodeex1.cpp \
           ode/linearode1storderex1.cpp \
           ode/nonlinearode1storderex1.cpp
