#!/usr/bin/perl

if(@ARGV==0) {
    print STDERR "This utility takes a gtf annotation (STDIN) and reformats it into a more compact, quickly readable file (STDOUT). Only exons are taken into account.\n";
}

parse_command_line(source  =>{description=>'the content of the source field', default=>'IPSA'},
		   element =>{description=>'the element', default=>'exon'},
		   features=>{description=>'list of features to output', default=>'transcript_id'},
		   lim     =>{description=>'line limit for debugging', default=>0});

foreach $item (split /\,/, $features) {
    ($f, $g) = split /\:/, $item;
    $g = 'uniq' unless($g);
    $feature_list{$f} = \&$g;
    print STDERR "[function $g($f)]\n"; 
}

print STDERR "[<stdin";
while($line=<STDIN>) {
    chomp $line;
    ($chr, $src, $elem, $beg, $end, $trash, $str, $trash, $attr) = split /\t/, $line;
    next unless($elem eq $element && $line=~/\w/);
    %attr =  get_attributes($attr);
    $tid = $attr{'transcript_id'};
    $CHR{$tid}{$chr} = $STR{$tid}{$str} = 1;
    push @{$exons{$tid}}, [$beg, $end];
    foreach $f(keys(%feature_list)) {
            $tr_feature{$f}{$tid} = $attr{$f};
    }
    $tr_feature{source}{$tid} = $src;
    $n++;
    last if($n>$lim && $lim>0);
}
print STDERR "]\n";

foreach $tid(keys(%exons)) {
    $N++;
    progressbar($N, 0+keys(%exons), "Processing ");
    if(keys(%{$CHR{$tid}})==1 && keys(%{$STR{$tid}})==1) {
	($chr, $str) = (keys(%{$CHR{$tid}}), keys(%{$STR{$tid}}));
	@exon_list = sort {$a->[0]<=>$b->[0]} @{$exons{$tid}};
	for($i=0,$j=0; $i<@exon_list; $i++) {
	    $key = join("\t", $chr, $source, 'exon',   $exon_list[$i]->[0],   $exon_list[$i]->[1], 0, $str, '.');
	    push @{$feature{$key}{transcript_id}}, $tid;
	    push @{$feature{$key}{position}}, relpos($i, @exon_list-1);
	    foreach $f(keys(%feature_list)) {
                push @{$feature{$key}{$f}}, $tr_feature{$f}{$tid} if($tr_feature{$f}{$tid});
            }
	    if($i>0) {
		$key = join("\t", $chr, $source, 'intron', $exon_list[$i-1]->[1], $exon_list[$i]->[0], 0, $str, '.');
		push @{$feature{$key}{transcript_id}}, $tid;   
		push @{$feature{$key}{position}}, relpos($i-1, @exon_list-2);
		foreach $f(keys(%feature_list)) {
                    push @{$feature{$key}{$f}}, $tr_feature{$f}{$tid} if($tr_feature{$f}{$tid});
            	}
	    }
	}
    }
    else {
	$trans_spliced++;
    }
}

print STDERR "[WARNING: $trans_spliced trans spliced transcripts excluded]" if($trans_spliced);

print STDERR "[>stdout";
foreach $key(sort keys(%feature)) {
    my %res=();
    foreach $f(sort keys(%feature_list)) {
        $res{$f} = $feature_list{$f}->(@{$feature{$key}{$f}}) if(@{$feature{$key}{$f}}>0);
    }
    print $key, "\t", set_attributes(%res), "\n";
}
print STDERR "]\n";

sub relpos {
    return('NA') unless(@_[0]>=0 && @_[1]>0);
    return(strand_c2i($str)<0 ? 1 - @_[0]/@_[1] : @_[0]/@_[1]);
}

sub abspos {
    return('I') if(@_[0]>0 && @_[0]<@_[1]);
    return(strand_c2i($str)>0 ? '5' : '3') if(@_[0]==0);
    return(strand_c2i($str)>0 ? '3' : '5');
}



################################################################################################################################################
# This routine parses the command line (@ARGV) according to the hash array specified in @_
# variable_name => {decription=>"blah", default=>"blah", ifunreadable=>"blah", ifabsent=>"blah", store=>"blah"}
# decription will be printed if @ARGV==0
# default is the default value of $variable_name
# ifunreadable will be printed if the file $variable_name is unreadable
# ifabsent will be printed if the $variable_name is empty
# store is TRUE/FALSE value of $variable_name
# array is a flag specifying that $variable_name is an variable (undef) array (array) or hash(hash)

sub parse_command_line {
    my %hash = @_;
    foreach $key(sort keys(%hash)) {
	my %param = %{$hash{$key}};
	($trash, $param{'default'}) = split /[=\n]/, `grep ^$key= makefile` if($param{'variable'});
	$$key = $param{'default'} if($param{'default'} ne undef);
	my $obligatory = $param{'ifunreadable'} || $param{'ifabsent'};
	if(@ARGV==0 && !$param{'variable'}) {
	    print STDERR "\t-$key", ($param{'store'} ? undef : " ..."), ", ", $param{'description'};
	    print STDERR ", default=$param{'default'}" if($param{'default'} ne undef);
	    print STDERR ", obligatory" if($obligatory);
	    print STDERR ", array=$param{'array'}" if($param{'array'});
	    print STDERR "\n";
	    next;
	}
	for(my $i=0;$i<@ARGV;$i++) {
	    if($ARGV[$i] eq "-$key") {
		if($param{'store'}) {
		    $$key = $param{'store'};
		}
		else {
		    if($param{'array'}) {
		    	push @{$key}, ($param{'array'} eq "hash" ? ($ARGV[++$i], $ARGV[++$i]) : $ARGV[++$i]);
		    }
		    else {
	    	        $$key = $ARGV[++$i];
		    }
		}
	    }
	}
	print STDERR "[WARNING: $key=$$key]\n" if($$key ne $param{'default'} && !$obligatory);
	die("ERROR: $key=$$key : $param{'ifunreadable'}\n") if($param{'ifunreadable'} && ! -r $$key);
	die("ERROR: $key=$$key : $param{'ifabsent'}\n") if($param{'ifabsent'} && $$key eq undef);
    }
    exit(1) if(@ARGV==0);
}


################################################################################################################################################
# get attribute field gtf style; input = string, output = hash
sub get_attributes {
    my %res = ();
    @_[0]="@_[0];" unless(@_[0]=~/;$/);
    while(@_[0]=~/([\w\_\d]+)\s*([\"\w\:\_\,\d\-\.]*)\;/g) {
	$res{$1} = $2;
	$res{$1} =~ s/^\"(.*)\"$/$1/;
    }
    return(%res);
}

# set attribute field gtf style input = hash, output = string
sub set_attributes {
    my %hash = @_;
    my @out=();
    foreach $key(sort keys(%hash)) {
        push @out, "$key \"$hash{$key}\";";
    }
    return(join(" ", @out));
}

# same but indexfile-like, not gtf-like
sub get_features {
    my %res = ();
    while(@_[0]=~/([\w\_\d]+)\=\"{0,1}(.*?)\"{0,1}\;/g) {
	$res{$1} = $2;
	#$res{$1} =~ s/\W//g;
    }
    return(%res);
}

sub set_features {
    my %hash = @_;
    my @out=();
    foreach $key(sort keys(%hash)) {
        push @out, "$key=\"$hash{$key}\";";
    }
    return(join(" ", @out));
}

################################################################################################################################################
# This routine generates and prints a command for the makefile as follows
#
# output : input depend (script if script_required) 
# 	script before input between output after
# endpoint :: output

sub update_path {
    my @res = ();
    foreach $item(@_) {
	if($item =~ /\.([A-Z]\d{2,2})\./) {
	    my @arr = split /\//, $item;
	    $arr[-1]="$1/$arr[-1]";
	    push @res, join("/", @arr);
	}
	else {
	    push @res, $item;
	}
    }
    return(@res);
}

sub makedir {
    my %dirs=();
    foreach my $filename(@_) {
        my @array = split /\//, $filename; 
	pop(@array);
        next unless(@array>0);
        $dirs{join("/", @array)}++;
    }
    print "\tmkdir -p ", join(" ", keys(%dirs)), "\n" if(keys(%dirs)>0);
}

sub make {
    my %param = @_;
    my %input  = %{$param{'input'}};
    my %output = %{$param{'output'}}; 
    push @{$param{'depend'}}, " $param{'script'}" if($param{'script_required'});
    $param{'script'} = "perl Perl/$param{'script'}" if($param{'script'}=~/\.pl$/);
    $param{'script'} = "Rscript R/$param{'script'}" if($param{'script'}=~/\.r$/);
    print join(" ", values(%output))," : ",join(" ", values(%input), @{$param{'depend'}}), "\n";
    makedir(values(%output));
    print "\ttouch ", join(" ", values(%output)),"\n" if($param{'touch'});
    print "\t$param{'script'} ",join(" ", $param{'before'}, %input, $param{'between'}, %output, $param{'after'})," \n";
    print "$param{'endpoint'} :: ", join(" ", values(%output)), "\n" if($param{'endpoint'});
    print "rm-$param{'endpoint'} ::\n\t rm -f ", join(" ", values(%output)), "\n" if($param{'endpoint'});
    print "touch ::\n\ttouch ", join(" ", values(%output)),"\n";
}

sub make2 {
    my %param = @_;
    my %inputs  = %{$param{'inputs'}};
    my %outputs = %{$param{'outputs'}};
    $param{'script'} = "perl Perl/$param{'script'}" if($param{'script'}=~/\.pl$/);
    foreach $key(keys(%{$param{'outputs'}})) { 
	print join(" ", values(%{$param{'outputs'}{$key}})),' '; 
    }
    print ":";
    foreach $key(keys(%inputs)) { 
	print join(" ", keys(%{$inputs{$key}})),' '; 
    }
    print "\n";
    foreach $key(keys(%outputs)) {
        makedir(values(%{$outputs{$key}}));
    }
    print "\t$param{'script'} $param{'before'} ";
    foreach $key(keys(%inputs)) {
	foreach $input(keys(%{$inputs{$key}})) {
	    print "$key $input $inputs{$key}{$input} ";
	}
    }
    print $param{'between'}, ' ';
    foreach $key(keys(%outputs)) {
	foreach $output(keys(%{$outputs{$key}})) {
	    print "$key $output $outputs{$key}{$output} ";
	}
    }
    print $param{'after'},"\n";
    return unless($param{'endpoint'});
    print "$param{'endpoint'} :: ";
    foreach $key(keys(%outputs)) {
	print join(" ", values(%{$outputs{$key}})), "\n";
    }
}

################################################################################################################################################
# returns formatted fraction (inclusion)/(inclusion + exclusion) if the denomintor > mincount; NA otherwise
sub frac {
    my ($INC, $EXC) = @_;
    my $TOT = $INC + $EXC;
    return($TOT > $mincount ?  sprintf("%.5f", $INC/$TOT) + 0 : "NA");
}
#
################################################################################################################################################
#
# strand char (+,.,-) to integer (1, 0 -1)
sub strand_c2i {
    return 1  if(@_[0] eq "+" || @_[0] eq "1");
    return -1 if(@_[0] eq "-" || @_[0] eq "-1");
    return 0;
}
#
#strand integer (1, 0 -1) to char (+,.,-)
sub strand_i2c {
    return("+") if(@_[0] eq 1);
    return("-") if(@_[0] eq -1);
    return '.';
}
#
################################################################################################################################################

sub progressbar {
    my ($current, $last, $message) = @_;
    $WCHAR =  80;
    my $width = 2**int(log($WCHAR)/log(2));
    my $i;
    if(int(($width*($current-1))/$last) < int(($width*$current)/$last)) {
        my $k = int(($width*$current)/$last);
        print STDERR "\r$message\[";
        for($i=0;$i<$k;$i++) {print STDERR "=";}
        print STDERR ">" if($k<$width);
        for($i++;$i<$width;$i++) { print STDERR " ";}
        print STDERR "\] ",($current/$last < 0.1 ? " " : ""), int(100*$current/$last),"%";
    }
    print STDERR "\n" if($current==$last);
}

################################################################################################################################################
sub sum {
    my $s=0;
    foreach my $x(@_) {$s+=$x;}
    return($s);
}

sub sumt {
    my $s=0;
    foreach my $x(@{@_[0]}) {
	$s++ if($x>@_[1]);
    }
    return($s);
}

sub avg {
    return("NA") unless(@_>0);
    return(sum(@_)/@_);
}

sub average {
    return(sprintf("%.2lf", avg(@_)));
}

sub min {
    return((sort{$a<=>$b}@_)[0])
}

sub max {
    return((sort{$a<=>$b}@_)[-1])
}

sub uniq {
    my %f = ();
    foreach $z(@_) {
        $f{$z}=1 if($z);
    }
    return(join(",", sort keys(%f)));
}

sub only {
    my @a = sort {length($a) <=> length($b)} @_;
    return($a[0]);
}

sub aggstat {
    my $s = 0;
    my $c = 0;
    my $l = 0;
    foreach $val(@_) {
        $s+=$val;
        $c+=1;
        $l+=$val*log($val);
    }
    my $h = sprintf("%.2f", (log($s) - $l/$s)/log(2));
    return($s, $c, $h>0 ? $h : 0);
}

sub read_annotation {
    my $N=0;
    print STDERR "[<@_[0]";
    open FILE, @_[0] || die ('Cannot read annotation');
	while($line=<FILE>) {
    	chomp $line;
    	($chr, $source, $feature, $beg, $end, $score, $str, $frame, $group) = split /\t/, $line;
    	if($feature eq "intron") {
    	    $SJ{$chr}{$beg}{$end}{$str}++;
    	    $SS{$chr}{$beg}{$str}++;
    	    $SS{$chr}{$end}{$str}++;
	}
	$N++;
    }
    close FILE;
    print STDERR ", $N annotated introns]\n";
}


sub annot_status {
    my ($chr, $beg, $end, $str) = @_;
    return($SJ{$chr}{$beg}{$end}{$str} ? 3 : ($SS{$chr}{$beg}{$str} && $SS{$chr}{$end}{$str} ? 2 : ($SS{$chr}{$beg}{$str} || $SS{$chr}{$end}{$str} ? 1 : 0)));
}
