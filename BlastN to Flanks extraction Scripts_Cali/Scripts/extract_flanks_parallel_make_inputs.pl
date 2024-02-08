# /usr/bin/env perl 

use warnings; 
use strict;
use POSIX qw/ceil/; 

# Number of parallel tasks to split the inputs into 
my $chunks = 224; # 28 X 8 = 224, for Broadwell nodes
 
my $dataset = 'IS110_complete';
my $filterName = 'Ident95_E0'; 

my $filtered_hits = "./Output/$dataset\_$filterName\.bacterial.blast.filtered";
my $list = "./Output/$dataset\_$filterName\.bacterial.blast-TEMPLIST";



  
#########################################################

# Hits per chunk:
my $hits = `wc -l < $filtered_hits`; 
chomp $hits; 
$hits--; # subtract for header
print "$filtered_hits has $hits hits\n"; 

my $hits_per_chunk = ceil($hits / $chunks); 
print "Printing $hits_per_chunk hits to each of $chunks files\n";  

 
open (H, $filtered_hits) || die "$! $filtered_hits\n";
chomp (my $header = <H>); 
 
my $chunk_num = 1; 
my $printed_hits = 0; 

my $cc = 0;  
  
while (my $line = <H> ) {
	chomp $line;
	$cc++; 
	my ($is_id, @rest) = split('\t', $line);
	my $out = "$filtered_hits\-TEMP-$chunk_num";
	
	if ( $printed_hits == 0 ) {
		open (O, ">$out") || die "$! $out\n"; 
		print O "$line\n"; 
		$printed_hits++; 
		#print "IF   : $printed_hits $chunk_num $cc\n"; 
		#print "Writing $line to $out\n";  
	}
	elsif ( $printed_hits < $hits_per_chunk ) {
		print O "$line\n";
		$printed_hits++; 
		#print "ELSIF: $printed_hits $chunk_num $cc\n"; 
		#print "Writing $line to $out\n";
	}
	else {
		#print "ELSE:  printed hits was $printed_hits will now be reset to 1. Chunk num was $chunk_num and will now be increased by 1\n";  
		close O;
		$chunk_num++; 
		$out = "$filtered_hits\-TEMP-$chunk_num";
		open (O, ">$out") || die "$! $out\n"; 
		print O "$line\n";
		$printed_hits = 1;
		#print "Writing $line to $out\n"; 		
	} 	
} close H; close O;  

# Create list of IS files sorted largest to smallest
#`wc -l $filtered_hits-TEMP-IS* | sort -rnk 1 | sed '1d' | awk '{print \$1","\$2}' > $list`;

open (L, ">$list") || die "$! write $list\n";
for (my $i = 1; $i <= $chunk_num; $i++) {
	print L "$filtered_hits\-TEMP-$i\n"; 
} close L; 
