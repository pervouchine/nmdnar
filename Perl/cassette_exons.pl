#!/usr/bin/perl

while($line=<STDIN>) {
    chomp $line;
    ($chr, $trash, $element, $beg, $end, $trash, $str, $frame, $attr) = split /\t/, $line;
    $data{$element}{$chr}{$str}{$beg}{$end}++;
    $atad{$element}{$chr}{$str}{$end}{$beg}++;
}

foreach $chr(keys(%{$data{exon}})) {
    foreach $str(keys(%{$data{exon}{$chr}})) {
	foreach $x(keys(%{$data{exon}{$chr}{$str}})) {
	    foreach $y(keys(%{$data{exon}{$chr}{$str}{$x}})) {
		foreach $a(keys(%{$atad{intron}{$chr}{$str}{$x}})) {
		    foreach $b(keys(%{$data{intron}{$chr}{$str}{$y}})) {
			next unless($data{intron}{$chr}{$str}{$a}{$b});
			$id = join("_", $chr, $x, $y, $str);
			print join("\t", $chr, $x, $y, $id, 1, $str), "\n";
		    }
		}
	    }
	}
    }
}



