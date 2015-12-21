#ifndef STABLE_H
#define STABLE_H

// Precompiled header
#if defined __cplusplus


#include <QtGlobal>

#ifdef Q_OS_WIN
#include <Windows.h>
#include "getopt.h"
#else
#include <getopt.h>
#endif

#ifndef NOGUI
#include <QApplication>
#include <QProcessEnvironment>
#include <QGLFormat>
#endif

// Add C++ includes here
#include <iostream>
#include <vector>
#include <utility>

typedef double number;

#include <configuration.h>
using Configuration = Sceneparser::Configuration<number>;
#include <linesegment.h>
using LineSegment = Sceneparser::LineSegment<number>;
#include <particle2.h>
using Particle2 = Particlegenerator::Particle2<number>;
using Particles = std::vector<Particle2*>;
#include <vector2.h>
using Vector2 = Sceneparser::Vector2<number>;
using Matrix2 = Sceneparser::Matrix2<number>;

using number_pair = std::pair<number, number>;

#endif

#endif // STABLE_H

