db_path = fullfile(pwd, "example.db");
db_schema = fullfile(pwd, "example.sql");

if exist(db_path, 'file')
    delete(db_path);
end

sqlcreate(db_path, db_schema);
