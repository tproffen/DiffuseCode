#!/usr/bin/perl

@f90files=glob("*.f90");
for ($i=0; $i<=$#f90files; $i++) {
   $useline=`grep USE $f90files[$i]`;
   @lines=split(/\n/,$useline);
   $out="";
   for ($j=0; $j<=$#lines; $j++) {
     if ($lines[$j]=~/^\s*USE/) {
       $lines[$j]=~s/^\s*USE\s*//;
       $lines[$j]=~s/\s*$//;
       unless (-r "$lines[$j].f90") {
         $lines[$j]="../../lib_f90/$lines[$j]";
       }
       unless ($out=~/$lines[$j]/) {
         $out.=$lines[$j];
         $out.=".mod ";
      }
    }
  }

  if ($out) { print "$f90files[$i]: $out\n"; }
}
