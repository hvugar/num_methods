TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
TARGET = main
#QMAKE_CXXFLAGS += -O2
#QMAKE_CXXFLAGS_RELEASE = -O3
#QMAKE_CFLAGS += -O3
#QMAKE_CXXFLAGS += -Wunused
#QMAKE_CXXFLAGS += -Werror
#QMAKE_CXXFLAGS -= -Wall
QMAKE_CXXFLAGS += -Wno-unused-variable
#QMAKE_CXXFLAGS += -O0
#QMAKE_CXXFLAGS += -ffloat-store
#QMAKE_CXXFLAGS += -ffast-math
QMAKE_CXXFLAGS += -D__USE_MINGW_ANSI_STDIO

INCLUDEPATH += ../minimum
LIBS += -L../bin -lminimum

INCLUDEPATH += ../border
LIBS += -L../bin -lborder

#INCLUDEPATH += ../hyperbolic
#LIBS += -L../bin -lhyperbolic

INCLUDEPATH += ../parabolic
LIBS += -L../bin -lparabolic

INCLUDEPATH += ../rnfunction
LIBS += -L../bin -lrnfunction

#INCLUDEPATH += ../imaging
#LIBS += -L../bin -limaging

DESTDIR += ../bin

SOURCES += main.cpp \
    problem1/iproblem1.cpp \
    problem1/problem1L2.cpp \
    problem1/problem1L3.cpp \
    problem1/art_problem1.cpp \
    problem1/problem1L1.cpp \
    high_order/singledifequ.cpp \
    high_order/systemdifequ.cpp \
    problem1/iloadedheatequation.cpp \
    problem1/loadedheatequation.cpp \
    problem1/ibackwardloadedheateqauation.cpp \
    problem4/problem4.cpp \
    ivp/cauchyproblemmex.cpp \
    load_sys/slodenlcs.cpp \
    problem5/problem5ex1.cpp \
    load_sys/slodenlcsm.cpp \
    problem5/zetta0.cpp \
    problem5/zetta1.cpp \
    problem5/zetta2.cpp

HEADERS += \
    problem1/iproblem1.h \
    problem1/problem1L2.h \
    problem1/problem1L3.h \
    problem1/art_problem1.h \
    problem1/problem1L1.h \
    high_order/singledifequ.h \
    high_order/systemdifequ.h \
    problem1/iloadedheatequation.h \
    problem1/loadedheatequation.h \
    problem1/ibackwardloadedheateqauation.h \
    problem4/problem4.h \
    ivp/cauchyproblemmex.h \
    load_sys/slodenlcs.h \
    problem5/problem5ex1.h \
    load_sys/slodenlcsm.h \
    problem5/zetta0.h \
    problem5/zetta1.h \
    problem5/zetta2.h
