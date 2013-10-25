TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

include(defaults.pri)

HEADERS += \
    src/diamondSquare/diamondSquare.h \
    src/mesher/mesher.h

SOURCES += \
    src/main.cpp \
    src/diamondSquare/diamondSquare.cpp \
    src/mesher/mesher.cpp
