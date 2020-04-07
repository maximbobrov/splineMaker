TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11
LIBS += -lopenGL32 -lGLU32
SOURCES += \
        main.cpp
QMAKE_CXXFLAGS_RELEASE += -O3 -ffast-math  -msse -std=c++11


QMAKE_LFLAGS += -O3 -ffast-math  -msse -std=c++11



LIBS += -lopenGL32 -lGLU32 -lm
LIBS += -L$$PWD/my_lib -lglut32
