function [db_exists, is_empty] = sql_exists(db_path)
% sql_is_empty  Check if db already exists and if it is empty
%   [db_exists, is_empty] = sql_exists(db_path)
%
% Arguments:
% - db_path: the db file path
%
% Output:
% - db_exists: true if the db exists
% - is_empty: true if the db is empty

if ~exist(db_path, 'file')
    db_exists = false;
    is_empty = true;
    return;
end

try
    db = sqlite(db_path, 'readonly');
    tables_tab = sqlread(db, "sqlite_master");
    num_tables = height(tables_tab);
    close(db);
    db_exists = true;
    is_empty = (num_tables == 0);
catch
    db_exists = false;
    is_empty = true;
end

end