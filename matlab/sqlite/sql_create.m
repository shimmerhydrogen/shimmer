function [] = sql_create(db_path, db_schema)
% sql_create  Create or fill a db from schema
%   [] = sql_create(db_path, db_schema) 
%
% Arguments:
% - db_path: the db file path
% - db_schema: the db file schema
%
% Note:
%   if the db does not exists the function creates it and fill it
%   if the db already exists and is empty the function fill it
%   if the db already exists and is not empty the function has no effect

[db_exists, is_empty] = sql_exists(db_path);

if ~is_empty
    return;
end

if ~db_exists
    conn = sqlite(db_path, 'create'); % Create the database
else
    conn = sqlite(db_path, 'connect'); % fill the database
end

sql_schema = fileread(db_schema); % Read schema from file

cmds = split(sql_schema, ';'); % Get all commands

for k=1:length(cmds)
    cmd = cmds{k};
    execute(conn, cmd); % Execute all commands
end

close(conn); % Close connection

end