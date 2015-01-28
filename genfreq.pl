#!/usr/bin/perl -w

use strict;
use List::Util qw[min max];
use Switch;
$|++; #unbuffer stout

#=== CONFIG START

#MISSION PARAMETERS
our $duration = 5; #mission duration $duration = 2, 5, 10 years
our $arm = 5; #armlength $arm = 1, 2, 3, 5, 5 million kilometers
#requires file name '"doppler_".$arm.".0.dat"' in folder '"/"$duration."year/"' with Doppler frequencies fd1, fd2, fd3 (one column per day) in Hz
our $links = 6; #number of $links = 4, 6

our $min = 7e6; #minimum frequency in Hz
our $max = 20e6; #maximum frequency in Hz

our $mode = 1; #1: no permutation of master laser. 2: do one run for each master laser position (6 runs). 3: find optimum maser laser position for each segment.
our $permutation = 123; #123 to 321 (only for $mode=1)
our $scheme = 1; #locking scheme $scheme = 1 (A), 2 (B), 3 (C)

#GENETIC PARAMETERS
our $population = 200; #population size
our $generations = 20; #number of generations
our $percentage = 50; #number of individuals allowed to procreate in percent
our $mutatingprob = 10; #mutation probability in percent
our $mutatingbits = 16; #least significant bits to mutate
our $recombination = 60; #probability that two offset frequencies are combined percent instead of one of both being used unchanged
our $alive = 1; #all offspring needs to survive for at least 1 day. 0: yes. 1: no.

#DEBUG
our $printdebug = 1; #write iterations for each generation. 0: no. 1: yes.
our $debugfolder = "debug"; #folder for debug files

#OUTPUT FILE DIRECTORY
my $mydate = &datum_zeit();
our $filedir = $mydate."_mission$duration-arm$arm"."_pop$population-gen$generations-per$percentage-int$recombination-mut$mutatingprob.$mutatingbits-min$min-max$max-scheme$scheme/";

#=== CONFIG END

my $exec;
$exec = `mkdir $filedir`;
$debugfolder = $filedir."debug";
$exec = `mkdir $debugfolder`;

our $filename = "";
our $switchmaster = 0;
our $run = 1;

switch ($mode) {
    case 1 {$switchmaster = 0; $run = 1;}
    case 2 {$switchmaster = 0; $run = 6;}
    case 3 {$switchmaster = 1; $run = 1;}
}

our $permutations;
if ($switchmaster == 1) {$permutations = 6;} else {$permutations = 1;}

our $A;
our $B;
our $C;
our @BN111; #beat-note between lasers of S/C1 on S/C1
our @BN121; #beat-note between lasers of S/C1 and S/C2 on S/C1
our @BN122; #beat-note between lasers of S/C1 and S/C2 on S/C2
our @BN131;
our @BN133;
our @BN222;
our @BN333;
our @BN232;
our @BN233;

our $end = 0;

#LOAD FILE
my $dopplerfile = $duration."year/doppler_".$arm.".0.dat";
my @FILE;
my $i = 0;
our @D; #Doppler shifts: ${$D[1]}[$day], ${$D[2]}[$day], ${$D[3]}[$day]

open (FILE, "$dopplerfile") or die "Cannot open $dopplerfile for reading: $!\n";
@FILE = <FILE>;
close FILE;

my $row;

foreach $row (@FILE) {
    
    if ($row !~ /^#/) {
        $end++;
        my @ROW;
        $row =~ s/\n//g;
        @ROW = split("\ +|\t+",$row);
        my $value;
        $i = 0;
        foreach $value (@ROW) {
            $i++;
            push(@{$D[$i]}, $value);
        }
    }
	
}

for (my $thisrun = 0; $thisrun<$run; $thisrun++) {

    if ($run > 1) {
        switch ($thisrun) {
            case 0 {$A = 1; $B = 2; $C = 3;}
            case 1 {$A = 3; $B = 2; $C = 1;}
            case 2 {$A = 2; $B = 3; $C = 1;}
            case 3 {$A = 1; $B = 3; $C = 2;}
            case 4 {$A = 3; $B = 1; $C = 2;}
            case 5 {$A = 2; $B = 1; $C = 3;}
        }
    }
    
    #output file name:
    switch ($mode) {
        case 1 {$filename = $filedir."$permutation";}
        case 2 {$filename = $filedir."$A$B$C";}
        case 3 {$filename = $filedir."optimum";}
    }
    
    #initiate output files
    my $OutF = $filename.".txt";
    open(OUTPUT, ">$OutF") or die "Can't open or create file $OutF: $!\n";
    print OUTPUT "#DAY #BN11\@1 #BN12\@2 #BN12\@1 #BN13\@3 #BN13\@1 #BN22\@2 #BN33\@3 #BN23\@2 #BN23\@3\n";
    close(OUTPUT);
    
    my $BrkF = $filename.".brk";
    open(BREAKS, ">$BrkF") or die "Can't open or create file $OutF: $!\n";
    print BREAKS "#Breaks\n";
    close(BREAKS);
    
    my $FrqF = $filename.".frq";
    open(FRQ, ">$FrqF") or die "Can't open or create file $FrqF: $!\n";
    print FRQ "#Frequency plan\n";
    close(FRQ);
    
    my $LOCK13 = 0; #lock between master laser and 2nd laser on S/C1
    my $LOCK21 = 0;
    my $LOCK31 = 0;
    my $LOCK23 = 0;
    my $LOCK32 = 0;
    
    my $start = 1; #start at day 1
    my $from = -$max;
    my $to = $max;
    
    print "Simulation for $end days.\n\n";
    
    open(FRQ, ">>$FrqF") or die "Can't open or create file $FrqF: $!\n";
    print FRQ "#Simulation for $end days.\n# min = $min;\n# max = $max;\n# from = $from;\n# to = $to;\n";
    print FRQ "#day LOCK13 LOCK21 LOCK31 LOCK23 LOCK32 days nice\n";
    close(FRQ);
    
    while ($start < $end) #try to keep beat-notes within frequency range for as long as possible. stop when end of file.
    {
        my @F = ();
        my @Fmax = ();
        my $dimensions = 5;
        my @POP;
        my $maxnice = 0;
        
        for (my $p = 0; $p<$permutations; $p++) {
            if ($run==1) {
                if ($switchmaster == 1) {
                    switch ($p) {
                        case 0 {$A = 1; $B = 2; $C = 3;}
                        case 1 {$A = 3; $B = 2; $C = 1;}
                        case 2 {$A = 2; $B = 3; $C = 1;}
                        case 3 {$A = 1; $B = 3; $C = 2;}
                        case 4 {$A = 3; $B = 1; $C = 2;}
                        case 5 {$A = 2; $B = 1; $C = 3;}
                    }
                    } else {
                    switch ($permutation) {
                        case 123 {$A = 1; $B = 2; $C = 3;}
                        case 321 {$A = 3; $B = 2; $C = 1;}
                        case 231 {$A = 2; $B = 3; $C = 1;}
                        case 132 {$A = 1; $B = 3; $C = 2;}
                        case 312 {$A = 3; $B = 1; $C = 2;}
                        case 213 {$A = 2; $B = 1; $C = 3;}
                    }
                }
            }
            
            @POP = ();
            
            #GENERATE GENERATION 0
            for (my $p = 0;$p<$population;$p++) {
                my $last = 0;
                my $nice = 0;
                my $warning = 0;
                my @sign = (-1,1); #positive or negative
                while ($last-$start < $alive) {
                    $F[0] = $sign[int(rand 2)]*(rand($max-$min)+$min)/1e6;
					$F[1] = $sign[int(rand 2)]*(rand($max-$min)+$min)/1e6;
					$F[2] = $sign[int(rand 2)]*(rand($max-$min)+$min)/1e6;
					$F[3] = $sign[int(rand 2)]*(rand($max-$min)+$min)/1e6;
					$F[4] = $sign[int(rand 2)]*(rand($max-$min)+$min)/1e6;
                    ($last,$nice) = beatnotes($start, $F[0], $F[1], $F[2], $F[3], $F[4]);
                    $warning++;
                    if (int($warning/100000) == $warning/100000) {print "WARNING: CANNOT FIND ANY COMBINATION IN $warning TRYS!\n"}
                }
                
                #nice in row 7: $last-$start-$nice. $last-$start is number of stable days, $nice is the mean value of the distance of the beat-notes from the center
                push (@POP, [$F[0], $F[1], $F[2], $F[3], $F[4], $last, $nice, $last-$start-$nice]);
                
                #debug:
                if ($printdebug) {
                    my $debugfile = "$debugfolder/day$start-$A$B$C.g0";
                    open(DEBUG, ">>$debugfile") or die "Can't open or create file $debugfile: $!\n";
                    print DEBUG $start." ".$F[0]." ".$F[1]." ".$F[2]." ".$F[3]." ".$F[4]." ".($last-$start)."\n";
                close(DEBUG); }
            }
            
            @POP = sort { $a->[7] <=> $b->[7] } @POP; #sort by $last-$start-$nice
            print "Day $start \@$A$B$C: stable for ".($POP[$population-1][5]-$start)." days, nice: $POP[$population-1][6] ($POP[$population-1][0], $POP[$population-1][1], $POP[$population-1][2], $POP[$population-1][3], $POP[$population-1][4])\n";
            
            #EVOLUTION
            my $selection = int($population*$percentage/100); #take the best of the best
            
            for (my $e = 0;$e<$generations;$e++) { #number of generations
                
                my @POPnew = ();
                
                #produce the next generation
                for (my $p = 0;$p<$population-1;$p++) {
                    
                    my $last = 0;
                    my $nice = 0;
                    my @Fnew = ();
                    
                    my $male = (int(rand($selection))+$population-$selection); #choose a male partner
                    
                    while ($last-$start < $alive) {
                        
                        my $female = (int(rand($selection))+$population-$selection); #choose a female partner
                        
                        my $negativ;
                        my $cut;
                        my $cut1;
                        my $cut2;
                        my $femalev;
                        my $malev;
                        my $interaction;
                        my $mutation;
                        my $mut;
                        
                        #PROCREATE
                        for (my $g = 0;$g<5;$g++) { #for all 5 offset lock frequencies
                            
                            $interaction = int(rand(100));
                            if ($interaction < $recombination) { #Interaktionsrate in %
                                
                                #use most significant bits of the father
                                if ($POP[$male][$g] < 0) {$negativ = -1} else {$negativ = 1}
                                $femalev = abs($POP[$female][$g]) * 100000 +0.5;
                                $malev = abs($POP[$male][$g]) * 100000 +0.5;
                                
                                #where to cut the gene
                                $cut = int(rand 17)+1;
                                $cut1 = 2**($cut)-1;
                                $cut2 = 0b1111111111111111111111111111111 << $cut;
                                
                                #cut male and female gene
                                $femalev = $femalev & $cut1;
                                $malev = $malev & $cut2;
                                
                                #offspring
                                $Fnew[$g] = $femalev | $malev;
                                $Fnew[$g] = $Fnew[$g] / 100000 * $negativ;
                                
							} else {
								$Fnew[$g] = $POP[$male][$g];
							}
                            
                            $mutation = int(rand(100));
                            if ($mutation < $mutatingprob) { #mutation rate in %
                                if ($Fnew[$g] < 0) {$negativ = -1} else {$negativ = 1}
                                $mut = 1 << int(rand $mutatingbits); #mutation: flip a bit somewhere
                                $Fnew[$g] = abs($Fnew[$g]) * 100000;
                                $Fnew[$g] = $Fnew[$g] ^ $mut;
                                $Fnew[$g] = $Fnew[$g] / 100000 * $negativ;
                            }                            
                        }
                        ($last,$nice) = beatnotes($start, $Fnew[0], $Fnew[1], $Fnew[2], $Fnew[3], $Fnew[4]);
                    }
                    
                    push (@POPnew, [$Fnew[0], $Fnew[1], $Fnew[2], $Fnew[3], $Fnew[4], $last, $nice, $last-$start-$nice]);
                    
                    #debug:
                    if ($printdebug) {
                        my $debugfile = "$debugfolder/day$start-$A$B$C.g".($e+1);
                        open(DEBUG, ">>$debugfile") or die "Can't open or create file $debugfile: $!\n";
                        print DEBUG $start." ".$Fnew[0]." ".$Fnew[1]." ".$Fnew[2]." ".$Fnew[3]." ".$Fnew[4]." ".($last-$start)."\n";
						close(DEBUG);
					}
                }
                #make sure that the strongest individual from the last generation is part of the new generation
                push (@POPnew, [$POP[$population-1][0], $POP[$population-1][1], $POP[$population-1][2], $POP[$population-1][3], $POP[$population-1][4], $POP[$population-1][5], $POP[$population-1][6], $POP[$population-1][7]]);
                
                if ($printdebug) {
                    my $debugfile = "$debugfolder/day$start-$A$B$C.g".($e+1);
                    open(DEBUG, ">>$debugfile") or die "Can't open or create file $debugfile: $!\n";
                    print DEBUG $start." ".$POP[$population-1][0] ." ".$POP[$population-1][1] ." ".$POP[$population-1][2] ." ".$POP[$population-1][3] ." ".$POP[$population-1][4] ." ".($POP[$population-1][5]-$start)."\n";
					close(DEBUG);
				}
                
                @POP = @POPnew; #replace generation
                
                @POP = sort { $a->[7] <=> $b->[7] } @POP; #sort by $last-$start-$nice
                print "Generation ".($e+1).": stable for ".($POP[$population-1][5]-$start)." days, nice: $POP[$population-1][6] ($POP[$population-1][0], $POP[$population-1][1], $POP[$population-1][2], $POP[$population-1][3], $POP[$population-1][4])\n";
                
            }
            
            if ($POP[$population-1][7] > $maxnice) {
                $maxnice = $POP[$population-1][7];
                $permutation = $A.$B.$C;
                $Fmax[0] = $POP[$population-1][0];
                $Fmax[1] = $POP[$population-1][1];
                $Fmax[2] = $POP[$population-1][2];
                $Fmax[3] = $POP[$population-1][3];
                $Fmax[4] = $POP[$population-1][4];
                $Fmax[5] = $POP[$population-1][5];
                $Fmax[6] = $POP[$population-1][6];
            }

        }


        switch ($permutation) {
            case 123 {$A = 1; $B = 2; $C = 3;}
            case 321 {$A = 3; $B = 2; $C = 1;}
            case 231 {$A = 2; $B = 3; $C = 1;}
            case 132 {$A = 1; $B = 3; $C = 2;}
            case 312 {$A = 3; $B = 1; $C = 2;}
            case 213 {$A = 2; $B = 1; $C = 3;}
        }

        
        #write the optimum frequency combination to file
        open(FRQ, ">>$FrqF") or die "Can't open or create file $FrqF: $!\n";
        print FRQ $start." ".$Fmax[0]." ".$Fmax[1]." ".$Fmax[2]." ".$Fmax[3]." ".$Fmax[4]." ".($Fmax[5]-$start)." ".$Fmax[6] ." ".$A.$B.$C."\n";
        close(FRQ);
        
        print "> ".($Fmax[5]-$start)." days \@$permutation!\n\n";
        
        #calculate the beat-notes
        &beatnotes($start, $Fmax[0], $Fmax[1], $Fmax[2], $Fmax[3], $Fmax[4]);
        
        #write the beat-notes to file
        &ausgabe($start,$Fmax[5]);
        
        $start = $Fmax[5]; #new starting point (day)
    }
    
}

sub beatnotes {
    
    my $start = shift;
    my $LOCK13 = shift;
    $LOCK13 = $LOCK13 * 1e6;
    my $LOCK21 = shift;
    $LOCK21 = $LOCK21 * 1e6;
    my $LOCK31 = shift;
    $LOCK31 = $LOCK31 * 1e6;
    my $LOCK23 = shift;
    $LOCK23 = $LOCK23 * 1e6;
    my $LOCK32 = shift;
    $LOCK32 = $LOCK32 * 1e6;
      
    my $i; #day to start at
    my $nice = 0;
    my $n = 0;
    for ($i=$start;$i<$end;$i++)
    {
        switch ($scheme) {
            case 1 {
                #LOCKING SCHEME A
                $BN111[$i] = $LOCK13;
                $BN122[$i] = $LOCK21;
                $BN121[$i] = 2*${$D[$A]}[$i]+$LOCK21;
                $BN133[$i] = $LOCK31;
                $BN131[$i] = 2*${$D[$C]}[$i]+$LOCK31;
                $BN222[$i] = $LOCK23;
                $BN333[$i] = $LOCK32;
                $BN232[$i] = ${$D[$A]}[$i]+$LOCK21+$LOCK23-($LOCK13+${$D[$C]}[$i]+$LOCK31+$LOCK32+${$D[$B]}[$i]);
                $BN233[$i] = ${$D[$A]}[$i]+$LOCK21+$LOCK23+${$D[$B]}[$i]-($LOCK13+${$D[$C]}[$i]+$LOCK31+$LOCK32);
                
                if ($links < 6) {
					$BN222[$i] = $min+0.5;
					$BN333[$i] = $min+0.5;
					$BN232[$i] = $min+0.5;
					$BN233[$i] = $min+0.5;
				}
            }
            case 2 {
                #LOCKING SCHEME B
                $BN111[$i] = ${$D[$A]}[$i]+${$D[$B]}[$i]+${$D[$C]}[$i]+$LOCK21+$LOCK23+$LOCK32+$LOCK31+$LOCK13;
                $BN122[$i] = $LOCK21;
                $BN121[$i] = 2*${$D[$A]}[$i]+$LOCK21;
                $BN133[$i] = 2*${$D[$C]}[$i]+$LOCK13;
                $BN131[$i] = $LOCK13;
                $BN222[$i] = $LOCK23;
                $BN333[$i] = $LOCK31;
                $BN232[$i] = 2*${$D[$B]}[$i]+$LOCK32;
                $BN233[$i] = $LOCK32;
                
            }
            case 3 {
				#LOCKING SCHEME C
                $BN111[$i] = $LOCK13;
                $BN121[$i] = 2*${$D[$A]}[$i]+$LOCK21;
                $BN131[$i] = ${$D[$A]}[$i]+${$D[$B]}[$i]+${$D[$C]}[$i]+$LOCK21+$LOCK23+$LOCK32+$LOCK31-$LOCK13;
                $BN222[$i] = $LOCK23;
                $BN122[$i] = $LOCK21;
                $BN232[$i] = 2*${$D[$B]}[$i]+$LOCK32;
                $BN333[$i] = $LOCK31;
                $BN133[$i] = $LOCK13+${$D[$C]}[$i]-${$D[$A]}[$i]-$LOCK21-$LOCK23-${$D[$B]}[$i]-$LOCK32-$LOCK31;
                $BN233[$i] = $LOCK32;
            }
        }
        
        $nice = $nice + (abs(abs($BN111[$i])-(($min+$max)/2))*abs(abs($BN111[$i])-(($min+$max)/2)) + abs(abs($BN122[$i])-(($min+$max)/2))*abs(abs($BN122[$i])-(($min+$max)/2)) + abs(abs($BN121[$i])-(($min+$max)/2))*abs(abs($BN121[$i])-(($min+$max)/2)) + abs(abs($BN133[$i])-(($min+$max)/2))*abs(abs($BN133[$i])-(($min+$max)/2)) + abs(abs($BN131[$i])-(($min+$max)/2))*abs(abs($BN131[$i])-(($min+$max)/2)) + abs(abs($BN222[$i])-(($min+$max)/2))*abs(abs($BN222[$i])-(($min+$max)/2)) + abs(abs($BN333[$i])-(($min+$max)/2))*abs(abs($BN333[$i])-(($min+$max)/2)) + abs(abs($BN232[$i])-(($min+$max)/2))*abs(abs($BN232[$i])-(($min+$max)/2)) + abs(abs($BN233[$i])-(($min+$max)/2))*abs(abs($BN233[$i])-(($min+$max)/2)))/9/(($max-$min)/2)/(($max-$min)/2); #mean distance from frequencies to center. border areas are more relevant. maximum is 1. smaller is better.
        $n++;
        
        if (abs($BN111[$i]) > $max || abs($BN111[$i]) < $min || abs($BN122[$i]) > $max || abs($BN122[$i]) < $min || abs($BN121[$i]) > $max || abs($BN121[$i]) < $min || abs($BN133[$i]) > $max || abs($BN133[$i]) < $min || abs($BN131[$i]) > $max || abs($BN131[$i]) < $min || abs($BN222[$i]) > $max || abs($BN222[$i]) < $min || abs($BN333[$i]) > $max || abs($BN333[$i]) < $min || abs($BN232[$i]) > $max || abs($BN232[$i]) < $min || abs($BN233[$i]) > $max || abs($BN233[$i]) < $min) {
            last;
        }
    }
        return($i,$nice/$n); #divide by number of days	

}


sub ausgabe {
    my $outstart = shift;
    my $outend = shift;
    my $OutF = $filename.".txt";
    open(OUTPUT, ">>$OutF") or die "Can't open or create file $OutF: $!\n";
    my $item;
    for ($i=$outstart;$i<$outend;$i++) {
        print OUTPUT "$i ".$BN111[$i]." ".$BN122[$i]." ".$BN121[$i]." ".$BN133[$i]." ".$BN131[$i]." ".$BN222[$i]." ".$BN333[$i]." ".$BN232[$i]." ".$BN233[$i]."\n";
    }
    close(OUTPUT);
    
    if ($i < $end) {
        my $BrkF = $filename.".brk";
        open(BREAKS, ">>$BrkF") or die "Can't open or create file $OutF: $!\n";
        print BREAKS "set arrow from $i.5,$min to $i.5,$max nohead linewidth 2 front\n";
        close(BREAKS);
    }
    
    return 1;
}


sub datum_zeit{
    
    my $p=$_[0];
    my %DATUM_ZEIT;
    my $timeparameter;
    my ($Sekunden, $Minuten, $Stunden, $Monatstag, $Monat, $Jahr, $Wochentag, $Jahrestag, $Sommerzeit) = localtime(time);
    $Monat+=1;
    $Jahrestag+=1;
    $Monat = $Monat < 10 ? $Monat = "0".$Monat : $Monat;
    $Monatstag = $Monatstag < 10 ? $Monatstag = "0".$Monatstag : $Monatstag;
    $Stunden = $Stunden < 10 ? $Stunden = "0".$Stunden : $Stunden;
    $Minuten = $Minuten < 10 ? $Minuten = "0".$Minuten : $Minuten;
    $Sekunden = $Sekunden < 10 ? $Sekunden = "0".$Sekunden : $Sekunden;
    $Jahr+=1900;
    
    $DATUM_ZEIT{'J'}=$Jahr;
    $DATUM_ZEIT{'M'}=$Monat;
    $DATUM_ZEIT{'D'}=$Monatstag;
    $DATUM_ZEIT{'h'}=$Stunden;
    $DATUM_ZEIT{'m'}=$Minuten;
    $DATUM_ZEIT{'s'}=$Sekunden;
    $DATUM_ZEIT{'T'}=$Wochentag;
    
    if ($p){
        $p=~s/(J|M|D|h|m|s|T)/$DATUM_ZEIT{$1}/g, while ($p=~/J|M|D|h|m|s|T/);
    }
    else{
        $p=$Jahr.$Monat.$Monatstag.'_'.$Stunden.$Minuten.$ Sekunden;
    }
    return $p;
}


$|++; #unbuffer (flush) STOUT

exit 0;
