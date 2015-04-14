TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

DEPENDPATH = ../minimum

HEADERS += ../minimum/methods.h \
           ../minimum/optimal.h \
           ../minimum/print.h
SOURCES += ../minimum/main.c \
           ../minimum/minimum.c \
           ../minimum/methods.c \
           ../minimum/optimal.c \
           ../minimum/optimal1.c \
           ../minimum/method_conj.c \
           ../minimum/method_grad.c \
           ../minimum/method_penalty.c \
           ../minimum/print.c \
           ../minimum/runga_kutta.c \
           ../minimum/sample_gradient.c \
           ../minimum/sample_functions.c \
           ../minimum/sample_penalty.c

