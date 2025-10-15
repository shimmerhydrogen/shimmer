package Shimmer;
use strict;
use warnings;

############################################################
use constant STATION_REMI   => 1;
use constant STATION_INJ    => 2;
use constant STATION_CONS   => 3;
use constant STATION_JUNC   => 4;
sub add_station {
    my ($db, $number, $type, $height) = @_;
    
    my $stat = $db->prepare(
        "INSERT INTO stations(s_number, s_name, t_type, s_height)
         VALUES (?, ?, ?, ?)"
    );

    my @record = ($number, "station_$number", $type, $height);
    $stat->execute(@record);
    $stat->finish;
    $db->commit;
}

############################################################
use constant PIPE_PLAIN     => 0;
sub add_pipe {
    my ($db, $from, $to) = @_;
    
    my $stat = $db->prepare(
        "INSERT INTO pipelines(p_name, s_from, s_to, p_type)
         VALUES (?, ?, ?, ?)"
    );

    my @record = ("pipe_${from}_${to}", $from, $to, PIPE_PLAIN);
    $stat->execute(@record);
    $stat->finish;
    $db->commit;
}

############################################################
sub add_remi_limit {
    my ($db, $station, %limit) = @_;
    
    my $stat = $db->prepare(
        "INSERT INTO limits_remi_wo VALUES (?, ?, ?, ?, ?)"
    );

    my @record = ($station, $limit{Lmin}, $limit{Lmax}, $limit{Pmin}, $limit{Pmax});
    $stat->execute(@record);
    $stat->finish;
    $db->commit;
}

############################################################
sub add_remi_const_profile {
    my ($db, $station, %profile) = @_;
    
    my $stat = $db->prepare(
        "INSERT INTO profiles_remi_wo VALUES (?, ?, ?)"
    );

    my @record = ($station, 0.0, $profile{Pset});
    $stat->execute(@record);
    $stat->finish;
    $db->commit;
}

############################################################
sub add_inj_limit {
    my ($db, $station, %limit) = @_;
    
    my $stat = $db->prepare(
        "INSERT INTO limits_injection_w(s_number, lim_Lmin, lim_Lmax,
            lim_Pmin, lim_Pmax) VALUES (?, ?, ?, ?, ?)"
    );

    my @record = ($station, $limit{Lmin}, $limit{Lmax}, $limit{Pmin}, $limit{Pmax});
    $stat->execute(@record);
    $stat->finish;
    $db->commit;
}

############################################################
sub add_inj_const_profile {
    my ($db, $station, %profile) = @_;
    
    my $stat = $db->prepare(
        "INSERT INTO profiles_injection_w VALUES (?, ?, ?, ?)"
    );

    my @record = ($station, 0.0, $profile{Pset}, $profile{Lset});
    $stat->execute(@record);
    $stat->finish;
    $db->commit;
}

############################################################
sub add_cons_limit {
    my ($db, $station, %limit) = @_;
    
    my $stat = $db->prepare(
        "INSERT INTO limits_conspoint_wo VALUES (?, ?, ?, ?, ?)"
    );

    my @record = ($station, $limit{Lmin}, $limit{Lmax}, $limit{Pmin}, $limit{Pmax});
    $stat->execute(@record);
    $stat->finish;
    $db->commit;
}

############################################################
sub add_cons_const_profile {
    my ($db, $station, %profile) = @_;
    
    my $stat = $db->prepare(
        "INSERT INTO profiles_conspoint_wo VALUES (?, ?, ?)"
    );

    my @record = ($station, 0.0, $profile{Lset});
    $stat->execute(@record);
    $stat->finish;
    $db->commit;
}

1;
