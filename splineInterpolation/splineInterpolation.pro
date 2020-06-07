TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11
LIBS += -lopenGL32 -lGLU32
SOURCES += \
        fileloader.cpp \
        globals.cpp \
        main.cpp \
        splineinterpolator.cpp \
        tools.cpp
#QMAKE_CXXFLAGS_RELEASE += -O3 -ffast-math  -msse -std=c++11
QMAKE_CXXFLAGS_RELEASE += -O2

#QMAKE_LFLAGS += -O3 -ffast-math  -msse -std=c++11
QMAKE_LFLAGS += -O2


LIBS += -lopenGL32 -lGLU32 -lm
LIBS += -L$$PWD/my_lib -lglut32

HEADERS += \
    fileloader.h \
    globals.h \
    splineinterpolator.h \
    tools.h
