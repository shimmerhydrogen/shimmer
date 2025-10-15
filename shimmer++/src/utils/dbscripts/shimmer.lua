if shimmer then return end

local posix = require 'posix';
local sql = require 'lsqlite3';

shimmer = {};

function shimmer.create_ndf(filename)
    local file = assert(io.open("../../../../sqlite/shimmer.sql", "r"));
    if (not file) then
        return nil;
    end

    if (posix.stat(filename)) then
        print("NDF already exists: " .. filename);
        return nil;
    end

    local handle = assert(io.popen("sqlite3 " .. filename, "w"));
    if (not handle) then
        return nil;
    end

    handle:write( file:read("*all") );

    handle:close();
    file:close();
end

function shimmer.open_ndf(filename)
    local ndf = sql.open(filename, sql.OPEN_READWRITE);
    return ndf;
end

function shimmer.close(ndf)
    sql.close(ndf);
end

shimmer.STATION_REMI = 1;
shimmer.STATION_INJ = 2;
shimmer.STATION_CONS = 3;
shimmer.STATION_JUNC = 4;

function shimmer.add_station(ndf, station, type, height)
    local stmt = ndf:prepare("INSERT INTO stations(s_number, \
        s_name, t_type, s_height) VALUES (?, ?, ?, ?)");
    
    stmt:bind_values(station, "station_" .. station, type, height);
    if (stmt:step() ~= sql.DONE) then
        print( ndf:errmsg() );
    end
    stmt:finalize();

    stmt = ndf:prepare("INSERT INTO gas_molar_fractions(s_number, \
        frac_CH4) VALUES (?, ?)");
    
    stmt:bind_values(station, 1.0);
    if (stmt:step() ~= sql.DONE) then
        print( ndf:errmsg() );
    end
    stmt:finalize();
end

function shimmer.add_remi_limit(ndf, station, limit)
    local stmt = ndf:prepare("INSERT INTO limits_remi_wo \
        VALUES (?, ?, ?, ?, ?)");
    
    stmt:bind_values(station, limit.Lmin, limit.Lmax, limit.Pmin, limit.Pmax);
    if (stmt:step() ~= sql.DONE) then
        print( ndf:errmsg() );
    end
    stmt:finalize();
end

function shimmer.add_remi_profile(ndf, station, profile)
    local stmt = ndf:prepare("INSERT INTO profiles_remi_wo \
    VALUES (?, ?, ?)");

    ndf:exec("BEGIN");
    for k, v in ipairs(profile) do
        stmt:bind_values(station, v[1], v[2]);
        if (stmt:step() ~= sql.DONE) then
            print( ndf:errmsg() );
        end
        stmt:reset();
    end
    stmt:finalize();
    ndf:exec("COMMIT");
end

function shimmer.add_inj_limit(ndf, station, limit)
    local stmt = ndf:prepare("INSERT INTO limits_injection_w(s_number, \
        lim_Lmin, lim_Lmax, lim_Pmin, lim_Pmax) VALUES (?, ?, ?, ?, ?)");
    
    stmt:bind_values(station, limit.Lmin, limit.Lmax, limit.Pmin, limit.Pmax);
    if (stmt:step() ~= sql.DONE) then
        print( ndf:errmsg() );
    end
    stmt:finalize();
end

function shimmer.add_inj_profile(ndf, station, profile)
    local stmt = ndf:prepare("INSERT INTO profiles_injection_w \
    VALUES (?, ?, ?, ?)");

    ndf:exec("BEGIN");
    for k, v in ipairs(profile) do
        stmt:bind_values(station, v[1], v[2], v[3]);
        if (stmt:step() ~= sql.DONE) then
            print( ndf:errmsg() );
        end
        stmt:reset();
    end
    stmt:finalize();
    ndf:exec("COMMIT");
end

function shimmer.add_cons_limit(ndf, station, limit)
    local stmt = ndf:prepare("INSERT INTO limits_conspoint_wo \
        VALUES (?, ?, ?, ?, ?)");
    
    stmt:bind_values(station, limit.Lmin, limit.Lmax, limit.Pmin, limit.Pmax);
    if (stmt:step() ~= sql.DONE) then
        print( ndf:errmsg() );
    end
    stmt:finalize();
end

function shimmer.add_cons_profile(ndf, station, profile)
    local stmt = ndf:prepare("INSERT INTO profiles_conspoint_wo \
    VALUES (?, ?, ?)");

    ndf:exec("BEGIN");
    for k, v in ipairs(profile) do
        stmt:bind_values(station, v[1], v[2]);
        if (stmt:step() ~= sql.DONE) then
            print( ndf:errmsg() );
        end
        stmt:reset();
    end
    stmt:finalize();
    ndf:exec("COMMIT");
end

function shimmer.set_sic(ndf, station, P, L)
    local stmt = ndf:prepare("INSERT INTO station_initial_conditions \
    VALUES (?, ?, ?)");
    stmt:bind_values(station, P, L);
    if (stmt:step() ~= sql.DONE) then
        print( ndf:errmsg() );
    end
    stmt:finalize();
end

function shimmer.add_pipe(ndf, from, to)
    local stmt = ndf:prepare("INSERT INTO pipelines(p_name, \
        s_from, s_to, p_type) VALUES (?, ?, ?, ?)");
    
    local pname = "pipe_" .. from .. "_" .. to;
    local PLAIN_PIPE = 0;
    stmt:bind_values(pname, from, to, PLAIN_PIPE);
    if (stmt:step() ~= sql.DONE) then
        print( ndf:errmsg() );
    end
    stmt:finalize();
end

function shimmer.set_pipe_params(ndf, from, to, diam, length, rough)
    local stmt = ndf:prepare("INSERT INTO pipe_parameters(p_name, s_from, \
        s_to, diameter, length, roughness) VALUES (?, ?, ?, ?, ?, ?)");

    local pname = "pipe_" .. from .. "_" .. to;
    stmt:bind_values(pname, from, to, diam, length, rough);
    if (stmt:step() ~= sql.DONE) then
        print( ndf:errmsg() );
    end
    stmt:finalize();
end

function shimmer.set_pic(ndf, from, to, G)
    local stmt = ndf:prepare("INSERT INTO pipe_initial_conditions \
    VALUES (?, ?, ?, ?)");
    local pname = "pipe_" .. from .. "_" .. to;
    stmt:bind_values(pname, from, to, G);
    if (stmt:step() ~= sql.DONE) then
        print( ndf:errmsg() );
    end
    stmt:finalize();
end

return shimmer;