
use strict;
#use warnings;                                                                                                                                                                                             

open ( IN, $ARGV[0] );

while ( my $line = <IN> ) {

    chomp $line;

    $new_name = $line;

    $new_name =~ s/counts.mod.txt.|.edgeR.DE_results//g;

    rename ( $line, $new_name );

}

