TEMPLATE  = lib
TARGET    = problem2P
CONFIG   -= app_bundle
#CONFIG  -= qt
CONFIG   += console
CONFIG   += shared
DESTDIR  = ../bin

DEFINES += PROBLEM2P_LIBRARY
#DEFINES += USE_IMAGING
#DEFINES += OLD_SOURCES

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

defined(USE_IMAGING) {
message("Using imaging library")
INCLUDEPATH += ../imaging
LIBS += -L../bin -limaging
QT       += core gui widgets
}

OBJECTS_DIR = release/.obj
MOC_DIR = release/.moc

HEADERS += problem2p_global.h \
           problem2p_common.h \
           problem2p_solver.h \
           problem2p_example.h

SOURCES += problem2p_common.cpp \
           problem2p_solver.cpp \
           problem2p_example.cpp

defined(OLD_SOURCES) {
    message("Using old sources")
#SOURCES += \
#    1d/problem2.cpp \
#    1d/iproblem2forward.cpp \
#    1d/iproblem2backward.cpp \
#    2d/iproblem2forward2d.cpp \
#    2d/iproblem2backward2d.cpp \
#    2d/cproblem2forward2d.cpp \
#    2d/cproblem2backward2d.cpp \
#    2d/iproblem2pibvp2d.cpp \
#    2d/abs/problem2forward2d.cpp \
#    2d/abs/problem2backward2d.cpp \
#    2d/dirakdelta.cpp \
#    2d/ex/expol.cpp \
#    2d/abs/ifunctional.cpp \
#    2d/abs/jfunctional.cpp \
#    2d/abs/pfunctional.cpp \
#    2d/ex/p2_article.cpp
#    #2d/abs/abstractproblem22d.cpp \
#    #2d/ex/problem22dex1.cpp \
#    #2d/ex/problem22dex2.cpp \
#    #2d/ex/problem22dex3.cpp \
#    #2d/ex/problem22dex4.cpp \
#    #2d/ex/problem22dex5.cpp \

#HEADERS += \
#    1d/problem2.h \
#    1d/iproblem2forward.h \
#    1d/iproblem2backward.h \
#    2d/iproblem2forward2d.h \
#    2d/iproblem2backward2d.h \
#    #2d/abs/abstractproblem22d.h \
#    #2d/ex/problem22dex1.h \
#    #2d/ex/problem22dex2.h \
#    #2d/ex/problem22dex3.h \
#    #2d/ex/problem22dex4.h \
#    #2d/ex/problem22dex5.h \
#    2d/cproblem2forward2d.h \
#    2d/cproblem2backward2d.h \
#    2d/iproblem2pibvp2d.h \
#    2d/abs/problem2forward2d.h \
#    2d/abs/problem2backward2d.h \
#    2d/dirakdelta.h \
#    2d/ex/expol.h \
#    2d/abs/ifunctional.h \
#    2d/abs/jfunctional.h \
#    2d/abs/pfunctional.h \
#    2d/ex/p2_article.h
}
