#include <QApplication>
#include <iostream>

#include "gui/CShimmerMainWindow.h"




int main(int argc, char **argv)
{
    QApplication app(argc, argv);

    CShimmerMainWindow w;
    w.setWindowTitle("Shimmer solver");
    w.show();

    return app.exec();
}