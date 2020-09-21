#! perl -w
##divide group based on 
#revision: 190702, replace all the sample id to fam id, specially for the samples that two id are diffrerent;


#use strict;
use Getopt::Long;

my ($prefix, $out, $fam, $dataset);

GetOptions  (
    "prefix=s"         => \$prefix,
    "fam=s"         => \$fam,
    "out=s"         => \$out,
    "dataset=s"       => \$dataset,
);

unless ($prefix && $fam) {
	die "Specify the prefix of input file and family file of sample list";
}

$out = $prefix unless ($out);

open(OUT1, "> $prefix.grp.stat") or die;
open(OUT2, "> $prefix.grp.det") or die;

##read group result file of dapc;
print scalar localtime;
print ":\tstart read fam file\n";
my %fam;
open(FAM, "$fam") or die;
while (<FAM>) {
    chomp;
    my @l=split /\s+/;
    $fam{$l[1]} = $l[0];
}
close FAM;

print ":\tstart read grp file\n";
my $grpfile = $prefix.".raw.grp";
open(GRP, "$grpfile") or die;
my %grp;
my %smpgrp;
my %grpsize;
my $grpnum;
my $validgrpnum;
while (<GRP>) {
    chomp;
    if (/"/) {
        $_ =~ s/"//g;
        my @l=split /\s+/;
        #$grp{$l[1]}{$l[0]} = 1;
        #$smpgrp{$l[0]} = $l[1];
        $grp{$l[1]}{$fam{$l[0]}} = 1; #replace the ID with FAM_ID
        $smpgrp{$fam{$l[0]}} = $l[1]; #replace the ID with FAM_ID
    }
    elsif (/\[\d+,\]/) {
        $_ =~ s/ +\[//g;
        $_ =~ s/\[//g;
        $_ =~ s/\]//g;
        $_ =~ s/\,//g;
        my @l=split /\s+/;
        $grpsize{$l[0]} = $l[1];
        $validgrpnum++ if ($l[1]>=54);
    }
    elsif (/size/) {my @l=split /\s+/; $grpnum = $l[1];}
}
close GRP;

print OUT1 "GRPNUM:\t$grpnum\nVALIDGRPNUM:\t$validgrpnum\n";


#print "GRPNUM:\t$grpnum\nVALIDGRPNUM:\t$validgrpnum\n";
##read fasta file;
print scalar localtime;
print ":\tstart read fasta file\n";
my $fastafile = $prefix.".fasta";
my %seq;
my $name;
my $len;
open(FASTA, "$fastafile") or die;
while (<FASTA>) {
    chomp;
    if (/>/) {$name = $_; $name =~ s/>//;}
    else {
        $len = length $_;
        my @base = split //, $_;
        for (1..$len) {
            $seq{$name}{$_} = $base[$_-1];
        }
    }
}
close FASTA;
##calc avr haplotype of group;
print scalar localtime;
print ":\tstart calc avr hap\n";
my %hap;
my %major;
for my $grp (keys %grp) {
    #next if ($grpsize{$grp}) < 54;
    my $tmp=$prefix.".$grp";
    open(TMP, "> $tmp") or die;
    for my $smp(keys %{$grp{$grp}}) {print TMP "$smp\n";}
    close TMP;
    
    system("plink --bfile $dataset --extract $prefix.list --keep-fam $tmp --out $tmp --freq > $tmp.log");
    my $n=0;
    open(FRQ, "$tmp.frq") or die;
    while (<FRQ>) {
        chomp;
        next if (/CHR/);
        $n++;
        my @l=split /\s+/;
        $major{$grp}{$n} = $l[4];
    }
    close FRQ;
    system("rm $tmp*");
}
##cmp sample with group major haplotype;
print scalar localtime;
print ":\tstart cmp smp with avr_hap\n";
my %cmp;
my $exceed;
my $totalvalid;
for my $grp (keys %grp) {
    next if ($grpsize{$grp} < 54);
    my @smp = keys %{$grp{$grp}};
    
    my ($dif, $valid, $dis);
    $cmp{$grp}{'exceed'} = 0;
    for my $smp (keys %{$grp{$grp}}) {
        for my $posi (sort {$a <=> $b} keys %{$seq{$smp}}) {
            next if ($seq{$smp}{$posi} ne '-');
            $valid++;
            $dif++ if ($seq{$smp}{$posi} ne $major{$grp}{$posi});
        }
        $dis = $dif / $valid if ($valid);
        $dis = 0 unless ($valid);
        $cmp{$grp}{'exceed'}++ if ($dis > 0.1);
    }
    $exceed += $cmp{$grp}{'exceed'};
    $totalvalid += $grpsize{$grp};
    print OUT2 "$grp\t$grpsize{$grp}\t$cmp{$grp}{'exceed'}\n";
}
my $exceedperc = $exceed / $totalvalid;
print OUT1 "TotalValid:Exceed:ExceedPerc\t$totalvalid:$exceed:$exceedperc\n";

#calc dis among groups
print scalar localtime;
print ":\tstart calc dis among avr_hap\n";
my %grpdis;
my %similar;
my %unsimilar;
my %grpdiscount;

my $similargrp = 0;
for my $grp (sort {$a <=> $b} keys %grp) {
    #next if ($grpsize{$grp} < 54);
    for my $grp2 (sort {$a <=> $b} keys %grp) {
        #next if ($grpsize{$grp2} < 54);
        my $dif = 0;
        if ($grp2 > $grp){
            for my $posi (keys %{$major{$grp}}) { $dif++ if ($major{$grp}{$posi} ne $major{$grp2}{$posi}); }
            $grpdis{$grp}{$grp2} = $dif / $len;
        }
        elsif ($grp eq $grp2)   {$grpdis{$grp}{$grp2} = "-";}
        else                    {$grpdis{$grp}{$grp2} = $grpdis{$grp2}{$grp};}
        
        if ($grpdis{$grp}{$grp2} eq '-') {}
        elsif ($grpdis{$grp}{$grp2} <= 0.1) {
            $grpdiscount{$grp}++;
            $similar{$grp}{$grp2} = 1;
            $similar{$grp2}{$grp} = 1;
        }
        elsif ($grpdis{$grp}{$grp2} > 0.1) {
            $unsimilar{$grp}{$grp2} = 1;
            $unsimilar{$grp2}{$grp} = 1;
        }
        
        #print OUT2 "$grpdis{$grp}{$grp2}\t";
        #printf OUT2 "%.2f\t", $grpdis{$grp}{$grp2};
        if ($grpdis{$grp}{$grp2} eq '-') {print OUT2 "-\t";}
        elsif ($grpdis{$grp}{$grp2} > 0.1) {printf OUT2 "%.2f\t", $grpdis{$grp}{$grp2};}
        elsif ($grpdis{$grp}{$grp2} <= 0.1) {
            printf OUT2 "%.2f\t", $grpdis{$grp}{$grp2};
            #push @{$similar{$grp}}, $grp2;
        }
    }
    print OUT2 "\n";
    $grpdiscount{$grp} = 0 unless ($grpdiscount{$grp});
    #$similargrp += $grpdiscount{$grp};
    
    my @similar = keys %{$similar{$grp}};
    my $similar = join ",", @similar;
    print OUT1 "GRP:SIZE:SIMILARGRP\t$grp:$grpsize{$grp}\t$grpdiscount{$grp}\t$similar\n";
}
##--------------merge similar groups
print scalar localtime;
print ":\tstart merge group\n";
my %merge;
my %merge2;
my %uncertain;
for my $grp (sort {$grpsize{$b} <=> $grpsize{$a}} keys %grp) {
    #next if ($grpsize{$grp} < 54 );
    next if ($merge{$grp});
    $merge{$grp} = $grp;
    $merge2{$grp}{$grp} = 1;
    
    my @grp = keys %{$similar{$grp}};
    unless (@grp) {next;}
    
    push @grp, $grp;
    
    for my $grp2 (keys %{$similar{$grp}}) {
        #next if ($grpsize{$grp2} < 54);
        next if ($merge{$grp2});
        my @grp2 = keys %{$similar{$grp2}};
        push @grp2, $grp2;
        
        my $tag=1;
        ###cmp similar collection of $grp and $grp2
        my %hash_a =  map {$_ => 1} @grp;
        #my %hash_b =  map {$_ => 1} @grp2;
        #my %merge_all = map {$_ => 1} @grp, @grp2;
        #my @a_only = grep {!$hash_b{$_}} @grp;
        my @b_only = grep {!$hash_a{$_}} @grp2;
        #my @common = grep {$hash_a{$_}} @grp2;
        #my @merge = keys (%merge_all);
        ###cmp dis of $grp/$grp2 and $grp2/other grp; judge grp assign
        if (@b_only >= 1) {
            for (@b_only) {$tag = 0 if ($grpdis{$grp2}{$_} < $grpdis{$grp2}{$grp});}
        }
        if ($tag) {$merge2{$grp}{$grp2} = 1; $merge{$grp2} = $grp;}
    }
}
######------------
print scalar localtime;
print ":\t output merge result\n";

print OUT1 "##merge result\n";
my %grpsize2;
for my $grp (sort {$grpsize{$b} <=> $grpsize{$a}} keys %merge2) {
    for my $grp2 (sort {$grpsize{$b} <=> $grpsize{$a}} keys %{$merge2{$grp}}) {$grpsize2{$grp} += $grpsize{$grp2};}
}

my $grpid = 0;
for my $grp (sort {$grpsize2{$b} <=> $grpsize2{$a}} keys %merge2) {
    $grpid++;
    print OUT1 "$grpid\t$grpsize2{$grp}\t";
    for my $grp2 (sort {$grpsize{$b} <=> $grpsize{$a}} keys %{$merge2{$grp}}) {$merge2{$grp2} = $grpid; print OUT1 "$grp2,";}
    print OUT1 "\n";
}

open(GRPOUT, "> $prefix.grp.new");

for my $smp (sort keys %seq) {
    print GRPOUT "$smp\t$merge2{$smpgrp{$smp}}\n";
}

#for my $grp (keys %grp) { for my $smp(keys %{$grp{$grp}}) {print GRPOUT "$smp\t$merge2{$grp}\n";} }

close OUT1;
close OUT2;
close GRPOUT;