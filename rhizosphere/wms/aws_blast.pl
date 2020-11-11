#!/usr/bin/env perl
#
#              INGLÊS/ENGLISH
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  http://www.gnu.org/copyleft/gpl.html
#
#
#             PORTUGUÊS/PORTUGUESE
#  Este programa é distribuído na expectativa de ser útil aos seus
#  usuários, porém NÃO TEM NENHUMA GARANTIA, EXPLÍCITAS OU IMPLÍCITAS,
#  COMERCIAIS OU DE ATENDIMENTO A UMA DETERMINADA FINALIDADE.  Consulte
#  a Licença Pública Geral GNU para maiores detalhes.
#  http://www.gnu.org/copyleft/gpl.html
#
#  Copyright (C) 2018  Universidade Estadual Paulista - UNESP
#
#  Universidade Estadual Paulista "Júlio de Mesquita Filho"
#  Laboratório de Bioinformática
#
#  Rafael Correia da Silva
#  rcs.biotec@gmail.com
#

=head1 NAME

=head1 SYNOPSIS

=head1 ABSTRACT

=head1 DESCRIPTION
    
    Arguments:

=head1 AUTHOR

Rafael Correia da Silva E<lt>rcs.biotec@gmail.comE<gt>

Copyright (c) 2018 Universidade Estadual Paulista - UNESP

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut

use strict;
use warnings;

#!/usr/bin/env perl
#
# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
# This software/database is a "United States Government Work" under the
# terms of the United States Copyright Act.  It was written as part of
# the author's official duties as a United States Government employee and
# thus cannot be copyrighted.  This software/database is freely available
# to the public for use. The National Library of Medicine and the U.S.
# Government have not placed any restriction on its use or reproduction.
#
# Although all reasonable efforts have been taken to ensure the accuracy
# and reliability of the software and data, the NLM and the U.S.
# Government do not and cannot warrant the performance or results that
# may be obtained by using this software or data. The NLM and the U.S.
# Government disclaim all warranties, express or implied, including
# warranties of performance, merchantability or fitness for any particular
# purpose.
#
# Please cite the author in any work or product based on this material.
#
# ===========================================================================
#
use strict;
use warnings;
use URI::Escape qw(uri_escape);
use LWP::UserAgent;
use LWP::Debug;
use HTTP::Request;
use Getopt::Long;
use Pod::Usage;
use File::Slurp qw(read_file);
use autodie;

my ($host, $program, $db, $query) = ("")x4;
my $format_type = undef;
my $evalue = undef;
my $filter = undef;
my $gap_costs = undef;
my $word_size = undef;
my $threshold = undef;
my $matrix = undef;
my $reward = undef;
my $penalty = undef;
my $num_descriptions = undef;
my $num_alignments = undef;
my $hitlist_size = undef;
my $job_title = undef;
my $cbs = undef;
my $num_threads = 1;
my $show_gis = 0;
my $verbose = 0;
my $dry_run = 0;
my $help_requested = 0;
GetOptions("host=s"         => \$host,
           "program=s"      => \$program,
           "db=s"           => \$db,
           "query=s"        => \$query,
           "format_type=s"  => \$format_type,
           "filter=s"       => \$filter,
           "evalue=s"       => \$evalue,
           "matrix=s"       => \$matrix,
           "gap_costs=s"    => \$gap_costs,
           "word_size=i"    => \$word_size,
           "threshold=s"    => \$threshold,
           "job_title=s"    => \$job_title,
           "reward=i"       => \$reward,
           "penalty=i"      => \$penalty,
           "cbs=i"          => \$cbs,
           "descriptions=i" => \$num_descriptions,
           "alignments=i"   => \$num_alignments,
           "hitlist_size=i" => \$hitlist_size,
           "num_threads=i"  => \$num_threads,
           "show_gis"       => \$show_gis,
           "verbose|v+"     => \$verbose,
           "dry_run"        => \$dry_run,
           "help|?"         => \$help_requested) || pod2usage(2);
pod2usage(-verbose=>2) if ($help_requested);
pod2usage(-message => "Missing hostname", -exitval=>1, -verbose=>2) unless (length $host);
pod2usage(-message => "Missing database", -exitval=>1, -verbose=>2) unless (length $db);
pod2usage(-message => "Missing program", -exitval=>1, -verbose=>2) unless (length $program);
pod2usage(-message => "Missing query", -exitval=>1, -verbose=>2) unless (length $query);
$verbose = 1 if ($dry_run and $verbose == 0);

$program = "blastn&MEGABLAST=on" if ($program =~ /megablast/i);

print "Submitting to $host: $program search of $db\n" if ($verbose >= 3);

my $encoded_query="";
my @contents = read_file($query);
foreach (@contents) {
    next if (/^$/);
    $encoded_query .= uri_escape($_);
}

# build the request
my $args = "CMD=Put&PROGRAM=$program&DATABASE=$db&QUERY=$encoded_query";
$args .= "&FORMAT_TYPE=$format_type" if (defined $format_type);
$args .= "&FILTER=$filter" if (defined $filter);
$args .= "&EXPECT=$evalue" if (defined $evalue);
$args .= "&MATRIX=$matrix" if (defined $matrix);
$args .= "&NUM_THREADS=$num_threads" if ($num_threads > 0);
$args .= "&GAPCOSTS=" . uri_escape($gap_costs) if (defined $gap_costs);
$args .= "&WORD_SIZE=$word_size" if (defined $word_size);
$args .= "&THRESHOLD=$threshold" if (defined $threshold);
$args .= "&NUCL_REWARD=$reward" if (defined $reward);
$args .= "&NUCL_PENALTY=$penalty" if (defined $penalty);
$args .= "&DESCRIPTIONS=$num_descriptions" if (defined $num_descriptions);
$args .= "&ALIGNMENTS=$num_alignments" if (defined $num_alignments);
$args .= "&HITLIST_SIZE=$hitlist_size" if (defined $hitlist_size);
$args .= "&COMPOSITION_BASED_STATISTICS=$cbs" if (defined $cbs);
$args .= "&JOB_TITLE=$job_title" if (defined $job_title);
$args .= "&NCBI_GI=T" if ($show_gis);
my $url = "http://$host/cgi-bin/blast.cgi";
print "$url\n" if ($verbose);

my $req = HTTP::Request->new(POST =>$url);
$req->content_type('application/x-www-form-urlencoded');
$req->content($args);
exit (0) if ($dry_run);

# get the response
my $ua = LWP::UserAgent->new(agent => "$0v1", keep_alive => 1);
$ua->default_header('Accept-Encoding' => scalar HTTP::Message::decodable());
$ua->show_progress(1) if ($verbose >= 3);
if ($verbose >= 2) {
    $ua->add_handler("request_send", sub { shift->dump; return; });
    $ua->add_handler("response_done", sub { shift->dump; return; });
}
my $response = $ua->request($req);

# parse out the request id
$response->decoded_content =~ /RID=(\w+)/;
my $rid=$1;
unless (defined $rid and length $rid) {
    my $err_msg = "Failed to obtain RID";
    $err_msg .= extract_error_msg($response->decoded_content);
    print STDERR "$err_msg\n";
    exit 3;
}

# Poll for results
$args = "CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";
while (1) {
    $req = HTTP::Request->new(GET=>"$url?$args");
    $response = $ua->request($req);
    print $response->decoded_content if $verbose;

    if ($response->decoded_content =~ /Status=ERROR/) {
        my $err_msg = "Search $rid failed";
        $err_msg .= extract_error_msg($response->decoded_content);
        print STDERR "$err_msg\n";
        exit 3;
    } elsif ($response->decoded_content =~ /Status=UNKNOWN/) {
        print STDERR "Search $rid expired.\n";
        exit 4;
    } elsif ($response->decoded_content =~ /Status=READY/) {
        last;
    }
    sleep(5);
}
$req = HTTP::Request->new(GET=>"$url?CMD=Get&RID=$rid");
$response = $ua->request($req);
print $response->decoded_content, "\n";
exit 0;

sub extract_error_msg
{
    my $content = shift;
    my $retval = "";
    if ($content =~ /class="error"/) {
        use List::MoreUtils qw(firstidx);
        my @content = split(/\n/, $content);
        my $idx = firstidx { /class="error"/ } @content;
        if ($idx != -1) {
            # Skip HTML mark up and extract error message
            for ($idx += 2;$content[$idx] !~ /<!--QBlastInfo/; $idx++) {
                $retval .= $content[$idx];
            }
        }
    }
    return $retval;
}


__END__

=head1 NAME

B<aws_blast.pl> - Submits BLAST searches via the URL API on the host specified

=head1 SYNOPSIS

aws_blast.pl -host <hostname> -program <program> -db <database> -query <query> [options]

    example: aws_blast.pl -host myhost -program blastp -db nr -query protein.fasta
    example: aws_blast.pl -host ec2-54-197-50-72.compute-1.amazonaws.com -program blastn -db pdbnt -query dna1.fasta
    example: aws_blast.pl -host ec2-54-197-50-72.compute-1.amazonaws.com -program megablast -db nt -query dna1.fasta

=head1 MANDATORY ARGUMENTS

=over

=item B<-host>

Host name to submit search to.

=item B<-program>

BLAST program to run, valid arguments include: blastn, megablast, blastp,
blastx, tblastn, tblastx.

=item B<-db>

BLAST database to search.

=item B<-query>

File containing FASTA sequence(s), GIs/accession numbers (one per line).

=back

=head1 OPTIONAL ARGUMENTS

=over

=item B<-format_type>

Format type to retrieve (valid arguments: html, xml, text, tabular, xml2, json2; default: html)

=item B<-filter>

Query filtering specification ("T" or "F", default: "T")

=item B<-evalue>

Expect value for BLAST search

=item B<-gap_costs>

A space separated pair of numbers representing the gap opening and extension
penalties

=item B<-word_size>

BLAST word size

=item B<-threshold>

BLAST word threshold

=item B<-cbs>

Composition based statistics value (0, 1, 2, or 3)

=item B<-job_title>

Title for BLAST search

=item B<-reward>

Reward for a match

=item B<-penalty>

Penalty for a mismatch

=item B<-descriptions>

Number of descriptions to show in the BLAST report

=item B<-alignments>

Number of alignments to show in the BLAST report

=item B<-hitlist_size>

Number of hits to keep

=item B<-show_gis>

Display NCBI GIs in BLAST results (default: false)

=item B<-num_threads>

Number of threads to run BLAST+ with (default: 1)

=item B<-verbose>, B<-v>

Produce verbose output (default: false)

=item B<-dry_run>

Do not run any commands, implies -verbose (default: false)

=item B<-help>, B<-h>, B<-?>

Displays this man page.

=back

=head1 EXIT CODES

    0 - success
    1 - missing command line arguments
    2 - invalid command line arguments
    3 - search failed
    4 - RID expired

=head1 AUTHOR

Christiam Camacho (camacho@ncbi.nlm.nih.gov)

=cut

