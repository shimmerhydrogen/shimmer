clear all;

db_path = fullfile(pwd, "example.db");
db_schema = fullfile(pwd, "example.sql");

%% reset test
if sql_exists(db_path)
    delete(db_path);
end

%% create db is empty
if ~sql_exists(db_path)
    sql_create(db_path, db_schema);
end

%% fill table
conn = sqlite(db_path, "connect");
users = { ...
     1, "io", "io_pwd", "io@mail.it"; ...
     2, "tu", "tu_pwd", "tu@mail.it"...
     };
users_tab = cell2table(users, ...
    "VariableNames", ["id", "username", "password", "email"] );
sqlwrite(conn, "users", users_tab);
close(conn);

%% read table
conn = sqlite(db_path, "readonly");
users_tab_1 = sqlread(conn, "users");
close(conn);

assert(height(users_tab_1) == 2);

%% reset database
[db_exists, is_empty] = sql_exists(db_path);

assert(db_exists);
assert(~is_empty);

if ~is_empty
    delete(db_path);
end
sql_create(db_path, db_schema);

%% read table
conn = sqlite(db_path, "readonly");
users_tab_2 = sqlread(conn, "users");
close(conn);

assert(height(users_tab_2) == 0);