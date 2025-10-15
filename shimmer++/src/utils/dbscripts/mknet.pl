#!/usr/bin/env perl
use strict;
use warnings;

use lib('.');
use Shimmer;

my $db_file = "test_gas_net.db";


############################################################
use DBI;
my $dbh = DBI->connect(
    "dbi:SQLite:dbname=$db_file", "", "", { 
        RaiseError => 1, AutoCommit => 0 }
    );

my %entry_limits = (
    Lmax => -300,
    Lmin => -10,
    Pmin => 6000000,
    Pmax => 8000000
);

my %exit_limits = (
    Lmin => 10,
    Lmax => 300,
    Pmin => 6000000,
    Pmax => 8000000
);


## Station 0    
Shimmer::add_station($dbh, 0, Shimmer::STATION_REMI, 100.0);
Shimmer::add_remi_limit($dbh, 0, %entry_limits);
Shimmer::add_remi_const_profile($dbh, 0, Pset => 7000000);

## Station 1
Shimmer::add_station($dbh, 1, Shimmer::STATION_JUNC, 200.0);

## Station 2
Shimmer::add_station($dbh, 2, Shimmer::STATION_CONS, 100.0);
Shimmer::add_cons_limit($dbh, 2, %exit_limits);
Shimmer::add_cons_const_profile($dbh, 2, Lset => 250);

## Pipes
Shimmer::add_pipe($dbh, 0, 1);
Shimmer::add_pipe($dbh, 1, 2);



$dbh->disconnect();


