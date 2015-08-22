TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.c \
    minimum.c \
    methods.c \
    method_grad.c \
    method_conj.c \
    method_prj_grad.c \
    print.c \
    sample_functions.c \
    sample_gradient.c \
    sample_penalty.c \
#    sample_grid.c \
    runga_kutta.c \
    method_penalty.c \
    method_grid.c \
    sample5.c \
    gradient.c

HEADERS += minimum.h \
    methods.h \
    print.h \
    method_grid.h \
    gradient.h \
    function.h
