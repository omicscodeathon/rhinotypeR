

<<com

List down the folder/dir which we wanted to loop through and store the basename which is the species
List down all fasta files within each directory
Extract the basename of the files excluding the directories
Read each file and capture the greater than sign (header), use cut and sed to extract the accession number
write the folder basename
com


folder=("RVA/" "RVB/" "RVC/")
for dir in ${folder[@]};
do
echo $dir
species=$(basename $dir)
echo $species
for file in $(ls $dir*.fasta)
do
genotype=$(basename $file | cut -f1 -d '.'|sed 's/genome//');echo $file;echo $genotype
accession=$(cat $file | grep '>' | cut -f1 -d ' ' | sed 's/>//')
echo $accession
out="$accession\t$genotype\t$species"
filename=Prototype.txt
echo -e $out >>$filename
done
done



