#!/bin/bash

IFS=$'\n'

#This scripts needs the Entrez Direct to work, see https://www.ncbi.nlm.nih.gov/books/NBK179288/ for download
#Otherwise, as much as possible is written in basic bash/perl in order to reduce the need for dependencies. This makes the code look a bit dirty, so feel free to adjust as desired.

#efetch path
edirectpath="/path/to/edirect"

#Path to the annotated kaiju output file
kaiju_output="/path/to/kaiju_summary_annotated.out"
#Which columns are taxonomy ID and species 
TaxidColumn="1"
SpeciesColumn="2"
SearchTerm="astrovir"


cat $kaiju_output | grep -i $SearchTerm | cut -f $TaxidColumn,$SpeciesColumn | perl -F'\t' -ane '$F[1] =~ s/[^0-9A-Za-z\n]/_/g; print "$F[0]\t$F[1]";' > virus_list.txt
cat $kaiju_output | grep -i $SearchTerm | cut -f $TaxidColumn,$SpeciesColumn > taxid_name.txt




for l in `cat virus_list.txt`
do
	name=`echo $l | cut -f 2`
	id=`echo $l | cut -f 1`
	F=$id"_"$name
	export F
	s="txid$id"
	echo $F
	mkdir $F
	cd $F
	$edirectpath/esearch -db nuccore -query "$s" | $edirectpath/efetch -format fasta > "$F"_fastas.fa
	
	cat "$F"_fastas.fa | perl -pe 's/^>(.*)\n/>$1HEADERLINE/g; s/\n//g; s/>/\n/g; s/HEADERLINE/\t/g;' | perl -pe 's/^\n//g;' > "$F"_fastas.txt
	cat "$F"_fastas.txt | perl -F'\t' -ane 'chomp $F[1]; $l=length($F[1]); @a=split(/ /, $F[0]); $acc=$a[0]; print "$l\t$acc\t$F[0]\n";' | sort -rnk 1 > "$F"_sequences.txt

	cd ..
done
