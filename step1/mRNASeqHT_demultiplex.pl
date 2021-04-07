#! /usr/bin/perl -w

#===============================================================================
#   Author: David Wang, @ Fluidigm inc
#
#   File: mRNASeqHT_demultiplex
#   Date: 01-09-2017
#   Version: 2.0.1
#
#   Usage:
#      mRNASeqHT_demultiplex.pl [options]
#
#      Try 'mRNASeqHT_demultiplex.pl -h' for more information.
#
#    Purpose: API of mRNASeqHT_demultiplex script will help you to demultiplex NGS data to individual single cells 
#			from Fluidigm's mRNASeq HT libs. the demultiplex are based on predefined row-barcodes in R1 reads.
#           The version does not require any non-core perl modules for processing.
#
#    Bugs: Please report bugs to technical support at Fluidigm
#			http://www.fluidigm.com 
#
#===============================================================================

use strict;


use Getopt::Std;
use vars qw/$opt_b $opt_i $opt_o $opt_h/;
use IO::File;

my $usage = "Usage: $0
	    -i input dir of fastq data
	    -o output dir of demultiplexed fastq data
	    [option]:
	    -h help
	    -b barcode file with tab-delimited format of <barcode_id> <barcode_seq>\n";


my %default_barcodes = (
	CACGTA	=>	'ROW01'	,
	CTCACA	=>	'ROW02'	,
	TGCATC	=>	'ROW03'	,
	TCAGAC	=>	'ROW04'	,
	CGATGT	=>	'ROW05'	,
	TACTGC	=>	'ROW06'	,
	ATGCTC	=>	'ROW07'	,
	CATCTG	=>	'ROW08'	,
	GACTCA	=>	'ROW09'	,
	AGATCG	=>	'ROW10'	,
	ATCAGC	=>	'ROW11'	,
	GCTACA	=>	'ROW12'	,
	CAGATC	=>	'ROW13'	,
	CACAGT	=>	'ROW14'	,
	TACGAG	=>	'ROW15'	,
	CGACTA	=>	'ROW16'	,
	GCATCT	=>	'ROW17'	,
	AGCACT	=>	'ROW18'	,
	ACACTG	=>	'ROW19'	,
	CGTAGA	=>	'ROW20'	,
	ACTCGA	=>	'ROW21'	,
	ACATGC	=>	'ROW22'	,
	CGTCAT	=>	'ROW23'	,
	TGTACG	=>	'ROW24'	,
	GCAGTA	=>	'ROW25'	,
	TCACGT	=>	'ROW26'	,
	ACGTCA	=>	'ROW27'	,
	CTCGAT	=>	'ROW28'	,
	ATCGTG	=>	'ROW29'	,
	GCTGAT	=>	'ROW30'	,
	GTCTAC	=>	'ROW31'	,
	CATGCT	=>	'ROW32'	,
	TAGCAC	=>	'ROW33'	,
	GTGCAT	=>	'ROW34'	,
	TAGTCG	=>	'ROW35'	,
	TCTCAG	=>	'ROW36'	,
	CTAGTC	=>	'ROW37'	,
	TCTAGC	=>	'ROW38'	,
	ATGACG	=>	'ROW39'	,
	GAGCTA	=>	'ROW40'	
);


my $cli = parse_cli();

if(defined $cli->{help}){
	die $usage,"\n";
}
my $fastq_data = get_fastq_data($cli->{fastq_dir});
my $barcodes = defined $cli->{barcode_file} ? get_barcode_data($cli->{barcode_file}) : \%default_barcodes;
my $out_dir = $cli->{out_dir};

system("mkdir $out_dir") if not -e $out_dir;

my %args = (
    barcode_len_in_R1	=> 6,
    fastq_data 			=> $fastq_data,
    barcodes			=> $barcodes,
    out_dir				=> $out_dir
);

demultiplex(\%args);

=head2 

	Title   : demultiplex
	Usage   : $data = demultiplex(args)
	Function: demultiplex from single-cell pool of column
	Returns : None
	Args    : args
    
=cut

sub demultiplex{
	my ($args) = @_;
	
	my $barcode_len_in_R1	= $args->{barcode_len_in_R1};
	my $fastq_data 			= $args->{fastq_data};
	my $barcodes 			= $args->{barcodes};
	my $out_dir 			= $args->{out_dir};
	
	my $revised_barcodes = get_barcode_hash($barcodes);
	
	my %summary = ();
	foreach my $sample(keys %$fastq_data){
		print  "Process  $sample ...\n";
		my $r1_fastq = $fastq_data->{$sample}->{R1};
		my $r2_fastq = $fastq_data->{$sample}->{R2};

		my $r1_in = IO::File->new("$r1_fastq") or die $!;
		my $r2_in = IO::File->new("$r2_fastq") or die $!;
		
		my %out_handles = ();
		
		foreach my $barcode(keys %$barcodes){
			my $id = $barcodes->{$barcode};
			my $cell_sample_name = join("_", $sample, $id);
			
			my $r1_sample_file = join("/", $out_dir, join("_", $cell_sample_name, "R1.fastq"));
			my $r2_sample_file = join("/", $out_dir, join("_", $cell_sample_name, "R2.fastq"));
			
			my $r1_sample_out = IO::File->new(">$r1_sample_file") or die $!;
			my $r2_sample_out = IO::File->new(">$r2_sample_file") or die $!;
			
			$out_handles{$cell_sample_name}{R1} = $r1_sample_out;
			$out_handles{$cell_sample_name}{R2} = $r2_sample_out;
			
		}
		
		
		my $undetermined_cell_sample_name =join("-",   $sample, "Undetermined");
		
		my $r1_undetermined_file = join("/", $out_dir, join("_", $undetermined_cell_sample_name, "R1.fastq"));
		my $r2_undetermined_file = join("/", $out_dir, join("_", $undetermined_cell_sample_name, "R2.fastq"));
			
		my $r1_undetermined_out = IO::File->new(">$r1_undetermined_file") or die $!;
		my $r2_undetermined_out = IO::File->new(">$r2_undetermined_file") or die $!;
			
		my $r1_barcode_out;
		my $r2_barcode_out;
		
		while (my $temp1_1 = <$r1_in>){
			
			my $temp1_2 = <$r1_in>;
    		my $temp1_3 = <$r1_in>;
    		my $temp1_4 = <$r1_in>;
    		
    		my $temp2_1 = <$r2_in>;
    		my $temp2_2 = <$r2_in>;
    		my $temp2_3 = <$r2_in>;
    		my $temp2_4 = <$r2_in>;
    	
    		my $bc = substr($temp1_2, 0, $barcode_len_in_R1);
    	
    		if (defined $revised_barcodes->{$bc}){
    			my $id = $revised_barcodes->{$bc};
    			my $cell_name = join("_", $sample, $id);
    			
    			$r1_barcode_out  = $out_handles{$cell_name}{R1};
    			$r2_barcode_out  = $out_handles{$cell_name}{R2};
    			
    			print $r1_barcode_out  join("", $temp1_1, $temp1_2,  $temp1_3, $temp1_4);	
    			print $r2_barcode_out  join("", $temp2_1, $temp2_2, $temp2_3, $temp2_4);
    			$summary{$sample}{$cell_name}++;
    		}
    		else{
    			print $r1_undetermined_out join("", $temp1_1, $temp1_2,  $temp1_3, $temp1_4);			
    			print $r2_undetermined_out join("", $temp2_1, $temp2_2, $temp2_3, $temp2_4);
    			$summary{$sample}{$undetermined_cell_sample_name}++;
    		}
	
    		
    	}
	}
	
	my $barcode_report = join("/", $out_dir, "demultiplex_report.xls");
	my $summary_out = IO::File->new(">$barcode_report") or die $!;	
	
	print $summary_out join("\t", "col_sample_name",  "cell_sample_name", "read_number"), "\n";
	
	foreach my $col_sample_name(sort keys %summary){
		foreach my $cell_sample_name(sort keys %{$summary{$col_sample_name}}){
			my $read_num = $summary{$col_sample_name}{$cell_sample_name};
			print $summary_out join("\t", $col_sample_name, $cell_sample_name, $read_num), "\n";
		}	
	}
	
}


=head2 

	Title   : get_barcode_data
	Usage   : $data = get_barcode_data(input_file)
	Function: read barcode data 
	Returns : hash of barcodes
	Args    : barcode file
    
=cut


sub get_barcode_data {
	my ($input) = @_;
	my $data = get_tab_delimited_data_from_file($input);
	my %barcodes = ();
	
	foreach my $row(@$data){
		my ($id, $barcode) = @$row;
		$barcodes{$barcode} = $id;
	}
	
	return(\%barcodes);
}

=head2 

	Title   : get_barcode_hash
	Usage   : $data = get_barcode_hash(barcodes)
	Function: revise barcodes to tolerate one seq error 
	Returns : hash of revised barcodes
	Args    : hash of orginal barcodes
    
=cut

sub get_barcode_hash{
	my ($barcodes) = @_;
    my %revised_barcodes = ();
	foreach my $barcode(keys %$barcodes){
	
		my $id = $barcodes->{$barcode};
		for(my $i =0; $i < length($barcode); $i++){
			foreach my $base (qw/A T G C N/){
				my $new_barcode = $barcode;
				substr($new_barcode, $i, 1) = $base;
				$revised_barcodes{$new_barcode} = $id;
			
			}
			
		}
	}
	
	return(\%revised_barcodes);
}


=head2 

	Title   : get_fastq_data
	Usage   : $data = get_fastq_data(input_dir)
	Function: get fastq data 
	Returns : hash of fastq files
	Args    : input dir
    
=cut

sub get_fastq_data {
	my ($dir) = @_;
    
    my %fastq = ();
    my %combined_fastq = ();
    
    opendir(my $dh, $dir) || die;
    while(my $file = readdir $dh) {
       next if $file !~ /fastq/;
       	my $command;
       
       if($file =~ /fastq\.gz/){
       	   my $file_str = join("/", $dir, $file);
       	  $command = "gunzip $file_str";
       	  `$command`;
       	  $file =~ s/\.gz//;
       }
       elsif($file =~ /fastq\.zip/ ){
       	 my $file_str = join("/", $dir, $file);
          $command = "unzip $file_str";
       	  `$command`;
       	  $file =~ s/\.zip//;
       }	
       
       my $sample_name;
       my $strand;
       
       if($file =~ /^(\S+)\_L00\d\_(R[1|2])/){
   	 		$sample_name = $1;
   	 		$strand = $2;
   		}
   		elsif($file =~ /^(\S+)\-L00\d\-(R[1|2])/){
   	 		$sample_name = $1;
   	 		$strand = $2;
   		}
   		elsif($file =~ /^(\S+)\_(R[1|2])/){
                $sample_name = $1;
                $strand = $2;
   		}
   		elsif($file =~ /^(\S+)\-(R[1|2])/){
                $sample_name = $1;
                $strand = $2;
   		}
   		elsif($file =~ /^(\S+)\_([1|2])\.fastq/){
                $sample_name = $1;
                $strand = $2;
   		}
   		elsif($file =~ /^(\S+)\-([1|2])\.fastq/){
                $sample_name = $1;
                $strand = $2;
   		}
   		else{
   			die "File name of $file is invalid\n ";
   		}
   		
   		$strand = join("", "R", $strand) if $strand !~ /R/;	
   		push @{$fastq{$sample_name}{$strand}}, join("/", $dir, $file);
   		
    }
    
    closedir $dh;
    
    
    foreach my $sample_name(keys %fastq){

    	if(not defined $fastq{$sample_name}{R1} or not defined $fastq{$sample_name}{R2}){
    		die "invalid fastq input, requre paired-end fastq data\n";
    	}
    	
    	my @r1_files = @{$fastq{$sample_name}{R1}};
    	my @r2_files = @{$fastq{$sample_name}{R2}};
    	
    	if(scalar @r1_files ==1 and scalar @r2_files ==1){
    		$combined_fastq{$sample_name}{R1} = $r1_files[0];
    		$combined_fastq{$sample_name}{R2} = $r2_files[0];
    	}
    	elsif(scalar @r1_files == scalar @r2_files){
    		my $combined_r1_file = join("/", $dir, join("_", $sample_name, "R1.fastq"));
    		my $combined_r2_file = join("/", $dir, join("_", $sample_name, "R2.fastq"));
    		my $command = join(" ", "cat " , join(" ", @r1_files),  "> $combined_r1_file"); 

    		system("$command");
    		$command =~ s/R1/R2/g;  
   	  	  	system("$command");
   	  	  
   	  	    $combined_fastq{$sample_name}{R1} =  $combined_r1_file;
    		$combined_fastq{$sample_name}{R2} =  $combined_r2_file;
   	  	  
    	}
    	elsif(scalar @r1_files != scalar @r2_files){
    		die "fastq file numbers in R1 and R2 of $sample_name are not equal\n";
    	}
    	
    }
	
	return(\%combined_fastq);
	
}


=head2 

	Title   : get_tab_delimited_data_from_file
	Usage   : $data = get_tab_delimited_data_from_file(input_file)
	Function: read tab-delimited data from file 
	Returns : array of tab-delimited data
	Args    : tab-delimited data
    
=cut

sub get_tab_delimited_data_from_file {
	my ($input) = @_;
	print "retrieve data from $input\n";
	my $in = IO::File->new($input) or die $!;
	my @data = ();
	
	while(<$in>){
		$_ =~ s/\s+$//;
		next if not $_;
		my @line = split("\t", $_);

		push @data, \@line;
		
	}

	return \@data;
}	


=head2 

	Title   : parse_cli
	Usage   : $data = parse_cli()
	Function: collect inputs
	Returns : hash of inputs
	Args    : None
    
=cut

sub parse_cli {
    getopts("i:b:o:h:");
    ($opt_i && $opt_o) ||
        die "Usage: $0
	    -i input dir of fastq data
	    -o output dir of demultiplexed fastq data
	    [option]:
	    -h help
	    -b barcode file with tab-delimited format of <barcode_id> <barcode_seq>\n";

    my %cli = (
                help					=> $opt_h,
		        fastq_dir 	            => $opt_i,
		        barcode_file			=> $opt_b,
		        out_dir					=> $opt_o);

    return \%cli;
}
