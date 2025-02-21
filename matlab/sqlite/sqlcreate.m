function [] = sqlcreate(db_path, db_schema)

if exist(db_path, 'file')
    return;
end

db = sqlite(db_path, 'create'); % Create or open the database
sql_schema = fileread(db_schema); % Read schema from file

cmds = split(sql_schema, ';');

for k=1:length(cmds)
    cmd = cmds{k};
    execute(db, cmd); % Execute all commands
end

close(db); % Close connection

end