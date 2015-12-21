#!/bin/bash --

wget http://www.libqglviewer.com/src/libQGLViewer-2.6.3.tar.gz
tar xzf libQGLViewer-2.6.3.tar.gz
rm libQGLViewer-2.6.3.tar.gz
cd libQGLViewer-2.6.3/QGLViewer/
qmake
make
sudo make install
cd ../../

