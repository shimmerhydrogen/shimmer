#include <iostream>
#include <string>

#include <QApplication>
#include <QtGui>
#include <QWidget>

#include "OpenXLSX.hpp"
#include "sol/sol.hpp"

using namespace OpenXLSX;

bool init_lua(sol::state& lua)
{
    lua.open_libraries(sol::lib::base);
    lua.create_table("config");
    sol::load_result script = lua.load_file("../share/shimmer_config.lua");
    if (!script.valid()) {
        std::cout << "Can't load config file" << std::endl;
        return false;
    }
    script();
    return true;
}

int main(int argc, char **argv)
{
    QApplication app(argc, argv);

    sol::state lua;
    if (not init_lua(lua))
        return 1;

    XLDocument doc;
    std::string path = lua["config"]["xlsx_input_file"];
    doc.open(path);
    std::string worksheet_name = lua["config"]["node_worksheet"];
    auto wks = doc.workbook().worksheet(worksheet_name);

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

