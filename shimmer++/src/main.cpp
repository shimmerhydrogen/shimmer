#include <iostream>
#include <string>

#include <QApplication>
#include <QtGui>
#include <QWidget>

#include "OpenXLSX.hpp"

using namespace OpenXLSX;

int main(int argc, char **argv)
{
    QApplication app(argc, argv);

    XLDocument doc;
    std::string path = "../../matlab/Mod_SHIMM_v01/INPUT_Pambour1.xlsx";
    doc.open(path);
    auto wks = doc.workbook().worksheet("PIPE");

    size_t nrows = 0;
    for (auto& row : wks.rows())
    {
        nrows++;
    }

    std::cout << nrows << std::endl;

    QWidget w;
    w.show();

    return app.exec();
}

