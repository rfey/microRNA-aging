#!/usr/bin/perl
use CGI;
use strict;

print "Content-type: text/html\n\n";
print CGI::start_html( -title => "Get 44 hr MicroRNA Expression Curve in Young and Old");

my $postData = getPostData();

my $mirBaseName = "";
if($postData->{mirBaseName}) {
    $mirBaseName = $postData->{mirBaseName};
}

# young rep 2 and 3, old rep 1 and 2: best-correlated replicates
my $dir = "/data/www/hendrixlab/html/data/circadian/microRNA_ATC/tables_RPMM";
my $youngRep2File = "/data/www/hendrixlab/html/data/circadian/microRNA_ATC/tables_RPMM/mirRPMMtable_youngRep2.txt";
my $youngRep3File = "/data/www/hendrixlab/html/data/circadian/microRNA_ATC/tables_RPMM/mirRPMMtable_youngRep3.txt";
my $oldRep1File = "/data/www/hendrixlab/html/data/circadian/microRNA_ATC/tables_RPMM/mirRPMMtable_oldRep1.txt";
my $oldRep2File = "/data/www/hendrixlab/html/data/circadian/microRNA_ATC/tables_RPMM/mirRPMMtable_oldRep2.txt";
my $gffFile = "/data/www/hendrixlab/html/data/circadian/microRNA_ATC/dme_dm6.gff3";
my $novelGffFile = "/data/www/hendrixlab/html/data/circadian/microRNA_ATC/dme_miRWoodsPredictions_youngAndOldNovel.gff";

my($known_mirNameToMirId) = readGeneList($gffFile);
my($novel_mirNameToMirId) = readGeneList($novelGffFile);
my($mirNameToMirId) = mergeMirIds($known_mirNameToMirId,$novel_mirNameToMirId); 
my @mirBaseNames = sort {lc $a cmp lc $b} keys(%{$mirNameToMirId});

print "<h1>44 hr expression plots for microRNAs</h1>";

print CGI::start_form(-method => "post",
		      -action => "/cgi-bin/getMirExpressionFigure_youngVsOldATC_RPMM.cgi");
print "<br>";

print "Select microRNA ID and click \"get figure\"\n<br>";
print "Many browsers will allow you to type the microRNA name once the menu is selected.<br>";
my $firstName = "Select MicroRNA Name";
print CGI::popup_menu(-name => "mirBaseName",
		      -values => \@mirBaseNames,
		      -default => $firstName);
print CGI::submit(-name => 'get figure',
		  -value => "get figure");
print CGI::end_form;
print "\n<br><hr>";

#####################
# BEGIN MAIN OUTPUT #
#####################

if($mirBaseName) {
    print "MicroRNA Expression Figure for $mirBaseName \n<br><br>";
    foreach my $mir (@{$mirNameToMirId->{$mirBaseName}}) {
	my($mirID,$mirName) = @{$mir};
	my $plotPrefix = $mirName . "_youngVsOld_replicatesPlot"; 
	my $outputFile = $dir . "/" . $plotPrefix . "_table.txt";
	my $youngVsOld_repsTable = parseFiles($youngRep2File,$youngRep3File,$oldRep1File,$oldRep2File,$outputFile,$mirID);
	# note the directory of the R Script must be readable from the web client
	#my $escapedGeneName = $mirName;
	#$escapedGeneName =~ s/\(/\\\(/g;
	#$escapedGeneName =~ s/\)/\\\)/g;
	#print "$mirName, $escapedGeneName\n";
	my $error = `Rscript /data/www/hendrixlab/html/data/circadian/microRNA_ATC/scripts/plotYoungVsOld_allReps_44hr_RPMM.R $youngVsOld_repsTable $mirName`;
	# note this line here must match with the output file of the plotYoungVsOld_allReps_noErrorBars_12hr.R script in /var/www/RScripts/
	print $error;
	my $imageName = "/data/circadian/microRNA_ATC/images_RPMM/" . $mirName . "_youngVsOld_replicates_44hrplot.pdf"; 
	my $imageFile = "/data/www/hendrixlab/html" . $imageName;
	print "<object data=\"$imageName\" type=\"application/pdf\" width=\"700\" height=\"740\"><a href=\"$imageName\">image file</a></object>";
    }
}

###############
# SUBROUTINES #
###############

sub readGeneList {
    my($gffFile) = @_;
    my %mirNameToMirId;
    my %baseDict;
    my %matureProducts;
    open(GFF,$gffFile) or die "Could not open $gffFile\n<br>";
    my %stored;
    while(<GFF>) {
	chomp;
	unless(/^#/){
	    my($chrom,$source,$type,$start,$stop,$score,$strand,$phase,$attribute) = split(/\t/);
	    if ($type eq "miRNA_primary_transcript"){
		my($baseID) = $attribute =~ /ID=(.*?);/;
		my($baseName) = $attribute =~ /Name=dme-(.*)/;
		$baseDict{$baseID} = $baseName;
	    }
	    if ($type eq "miRNA") {
		my($derivedFrom) = $attribute =~ /Derives_from=(.*)/;
		my($mirID) = $attribute =~ /ID=(.*?);/;
		my($mirName) = $attribute =~ /Name=(.*?);/;
		my($shortName) = $mirName =~ /dme-(.*)/;
		if ($mirID =~ /MIMAT.*/) {
		    $mirID = "${mirName}_${mirID}"; 
		}
		push(@{$matureProducts{$derivedFrom}},[$mirID,$shortName]);
	    }
	}
    }
    foreach my $key ( keys %baseDict ) {
	my $newkey = $baseDict{$key};
	my $val = $matureProducts{$key};
	$mirNameToMirId{$newkey} = $val;
    }
    return \%mirNameToMirId;
    close(GFF);
}

sub mergeMirIds {
    my($known_mirNameToMirId,$novel_mirNameToMirId) = @_;
    foreach my $key ( keys %{$novel_mirNameToMirId} ) {
        $known_mirNameToMirId->{$key} = $novel_mirNameToMirId->{$key};
    }
    return $known_mirNameToMirId;
}

sub getPostData {
    my %postData;
    foreach my $name ( CGI::param ) {
	$postData{$name} = CGI::param($name);
    }
    return \%postData;
}

sub readRepFile {
    my($mirName,$repFile) = @_;
    #print "$mirName\n$repFile\n<br>";
    my(@replicate);
    open(REP,$repFile) or die "Could not open $repFile\n";
    while(<REP>) {
	chomp;
	unless(/^#mirID/) {
	    my($mirId,$ZT0,$ZT4,$ZT8,$ZT12,$ZT16,$ZT20) = split(/\s+/);
	    #my($mirBaseID,$mirAcc) = split(/_/,$mirId);
	    if($mirName eq $mirId) {
		@replicate = ($ZT0,$ZT4,$ZT8,$ZT12,$ZT16,$ZT20);
		#print "$mirId,$ZT0,$ZT4,$ZT8,$ZT12,$ZT16,$ZT20\n";
		last;
	    }
	}
    }
    close(REP);
    return @replicate;
}

sub parseFiles {
    my($youngRep2File,$youngRep3File,$oldRep1File,$oldRep2File,$outputFile,$mirName) = @_;
    my @young2 = readRepFile($mirName,$youngRep2File);
    my @young3 = readRepFile($mirName,$youngRep3File);
    my @old1 = readRepFile($mirName,$oldRep1File);
    my @old2 = readRepFile($mirName,$oldRep2File);
    open(TABLE,">$outputFile") or die "Could not open $outputFile for writing.\n";
    my @ZT = (0,4,8,12,16,20);
    for(my $i=0; $i<@ZT;$i++) {
	print TABLE "$ZT[$i]\t$young2[$i]\t$young3[$i]\t$old1[$i]\t$old2[$i]\n"; 
    }
    close(TABLE);
    return $outputFile;
}
