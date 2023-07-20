/* This code is part of the SHIMMER project
 *
 * Politecnico di Torino, Dipartimento di Matematica (DISMA)
 * 
 * The authors (C) 2023 
 */

/* Note: these import functions exist only for convenience. Xlsx is an horrible
 * format to work with, a more robust input format will be needed. */

#include "OpenXLSX.hpp"
#include "infrastructure_graph.h"
#include "xlsx_io.h"

using namespace OpenXLSX;

#define EDGE_ROW_LENGTH 11
#define EDGE_BRANCH     0
#define EDGE_FROM       1
#define EDGE_TO         2
#define EDGE_LENGTH     3
#define EDGE_DIAMETER   4
#define EDGE_GRID_PTS   10

static int
vertices_from_xlsx(const XLWorksheet& wks, std::vector<vertex_properties>& vps)
{
    std::vector<XLCellValue> row_values;
    int parsed_rows = 0;
    for (auto& row : wks.rows())
    {

    }

    return 0;
}

static int
edges_from_xlsx(const XLWorksheet& wks, std::vector<edge_properties>& eps)
{
    std::vector<XLCellValue> row_values;
    int parsed_rows = 0;
    for (auto& row : wks.rows())
    {
        if (parsed_rows == 0) {
            parsed_rows = 1;
            continue; // skip header
        }

        row_values = row.values();

        if (row_values.size() < EDGE_ROW_LENGTH)
            break;

        if (row_values[EDGE_BRANCH].type() == XLValueType::Empty)
            break;
    
        edge_properties ep;
        ep.branch_num   = row_values[EDGE_BRANCH];
        ep.from         = row_values[EDGE_FROM];
        ep.to           = row_values[EDGE_TO];
        //if a cell contains a round number it gets interpreted as an integer
        //and the type conversion fails crashing the program.
        //ep.length       = row_values[EDGE_LENGTH];
        //ep.diameter     = row_values[EDGE_DIAMETER];
        //ep.grid_pts     = row_values[EDGE_GRID_PTS];
        ep.type         = edge_type::pipe; /* XXX: handle the other cases */

        eps.push_back( std::move(ep) );
    }

    return parsed_rows;
}

int
import_infrastructure_from_xlsx(sol::state& lua, infrastructure_graph& ig)
{
    XLDocument xlsx;

    auto path = lua["config"]["xlsx_input_file"];
    if (not path.valid()) {
        std::cerr << "Please set 'config.xls_input_file' in the configuration." << std::endl;
        return -1;
    }

    auto node_sheet_name = lua["config"]["node_worksheet"];
    if (not node_sheet_name.valid()) {
        std::cerr << "Please set 'config.node_worksheet' in the configuration." << std::endl;
        return -1;
    }

    auto edge_sheet_name = lua["config"]["edge_worksheet"];
    if (not edge_sheet_name.valid()) {
        std::cerr << "Please set 'config.edge_worksheet' in the configuration." << std::endl;
        return -1;
    }

    xlsx.open(std::string(path));
    auto node_sheet = xlsx.workbook().worksheet(std::string(node_sheet_name));
    auto edge_sheet = xlsx.workbook().worksheet(std::string(edge_sheet_name));

    std::vector<vertex_properties> vps;
    int err = vertices_from_xlsx(node_sheet, vps);

    std::vector<edge_properties> eps;
    err = edges_from_xlsx(edge_sheet, eps);

    std::cout << eps.size() << std::endl;

    return 0;
}