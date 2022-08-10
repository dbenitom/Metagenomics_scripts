# Original script by Viola Krukenberg (https://orcid.org/0000-0001-8369-8114)
# With modifications by David Benito Merino (https://orcid.org/my-orcid?orcid=0000-0002-1453-7330, https://github.com/dbenitom)

# This script detects CxxCH (heme-binding) motifs in predicted open reading frames from microbial genomes or metagenomes
# and imports the results in an anvi'o contigs database (https://anvio.org/).

# Path to anvio contigs.db:
CONTIGS_DB=
# Annotated genome name:
GENOME_ANNOTATION=
# Threads:
THREADS=10
# output
WORK_DIR=
#motif to search for
MOTIF="C[A-Z][A-Z]CH"
#short name of the motif
MOTIF_NAME="CXXCH"


####
# Command to extend aliases from bash profile
shopt -s expand_aliases

echo "## Step 1. Generate the fasta file with gene calls for the genome/metagenome:"
date
echo ""

source activate anvio-master; source ~/virtual-envs/anvio-master/bin/activate

anvi-get-sequences-for-gene-calls -c $CONTIGS_DB --get-aa-sequences -o ${WORK_DIR}/${GENOME_ANNOTATION}.faa

deactivate; source deactivate

echo "## Done"
date
echo ""

echo "## Step 2. Search motif in the gene calls:"
date
echo ""
cd $WORK_DIR
mkdir ${MOTIF_NAME}_scan

cp ${WORK_DIR}/${GENOME_ANNOTATION}.faa ${MOTIF_NAME}_scan

cd ${MOTIF_NAME}_scan

#modify annotated .faa file: each sequence in a single line
awk  '{ if ($1 ~ /^>/) { printf("\n%s\n", $0) } else { printf("%s", $0) }}' ${GENOME_ANNOTATION}.faa > modified_${GENOME_ANNOTATION}.faa

#find sequences containing the motif and output them to new file
grep -B 1 ${MOTIF} modified_${GENOME_ANNOTATION}.faa > ${MOTIF_NAME}_motif_${GENOME_ANNOTATION}.faa

#make a list with header (including annotation) of the sequences containing the motif
grep '^>' ${MOTIF_NAME}_motif_${GENOME_ANNOTATION}.faa | sed 's/^>//' > motifList_${GENOME_ANNOTATION}.txt

#make a list of only IDs of the motif containing sequences
grep '^>'  ${MOTIF_NAME}_motif_${GENOME_ANNOTATION}.faa | sed 's/>//' | sed 's/ .*//' > motifIds_${GENOME_ANNOTATION}.txt

#this extracts sequences using a list of IDs; script from Brandon; make sure script is in path
grep -A 1 -f motifIds_${GENOME_ANNOTATION}.txt modified_${GENOME_ANNOTATION}.faa > Sequences_${MOTIF_NAME}motif_${GENOME_ANNOTATION}.fasta

#make directory for count data
mkdir motifCountTables_${MOTIF_NAME}

#copy the file containing the sequences to the new directory; csplit will create many which should be stored in the directory for count data; could also move files here after cspli
cp Sequences_${MOTIF_NAME}motif_${GENOME_ANNOTATION}.fasta motifCountTables_${MOTIF_NAME}

cd  motifCountTables_${MOTIF_NAME}

#split the multifasta file into separate files each containing one sequence; did this to then count in the next step the motifs per sequence
csplit -f ${MOTIF_NAME} Sequences_${MOTIF_NAME}motif_${GENOME_ANNOTATION}.fasta '/>/' '{*}'

#for each fasta file from csplit output make a file with the ORF name, count the number of motifs
for File in ${MOTIF_NAME}*
  do
    echo "$File"
    grep '>' "$File" | sed 's/^>//' > "$File"_name

    #grep every motif (-o for only matching patterns) and make a list with all motifs founds
    grep -o ${MOTIF} "$File" >  "$File"_list

    #count the number of lines in the list of motifs found to know the number of motifs in the sequence
    wc -l "$File"_list > "$File"_count

    #make a table with the number of motifs and the ORF name, separate the columns by a tab
    paste "$File"_count  "$File"_name   | column -c 2 -s $'\t' > CountTable_"$File".tab
done


#concatenate the tables created for each file into one big table
cat CountTable*.tab > CountTableAll_${MOTIF_NAME}.tab

#use tail to remove first line (first file contains no data), use sed to remove file name after count (separated by a space) and add a header to the columns of the table
tail -n +2 CountTableAll_${MOTIF_NAME}.tab | sed 's/ .*\t/\t/' | sed '1i\motif#\tORF' > FinalCountTableAll_${MOTIF_NAME}.tab

#clean the directory by moving the in previous steps generated files to separte folders; could also delete them
mkdir Names
mv *_name Names/

mkdir Lists
mv *_list Lists/

mkdir Counts
mv *_count Counts/

mkdir CountTables
mv CountTable_* CountTables/

mkdir SequFiles
mv ${MOTIF_NAME}* SequFiles/

echo "## Done"
echo ""


echo "## Step 3. Parse and import to contigs database"
date
echo ""
# Column 1: gene calls
cut -f2 FinalCountTableAll_${MOTIF_NAME}.tab | tail --lines=+2 > ${WORK_DIR}/gene_calls.tmp
# Column 2: source
yes "CXXCH_scan" | head -n `cut -f2 ${WORK_DIR}/gene_calls.tmp | wc -l` > ${WORK_DIR}/source.tmp
# Column 3: accession
yes "CXXCH_count" | head -n `cut -f2 ${WORK_DIR}/gene_calls.tmp | wc -l` > ${WORK_DIR}/accession.tmp
# Column 4: function
cut -f1 FinalCountTableAll_${MOTIF_NAME}.tab | tail --lines=+2 > ${WORK_DIR}/function.tmp
# Column 5: e-value
yes "0" | head -n `cut -f2 ${WORK_DIR}/gene_calls.tmp | wc -l` > ${WORK_DIR}/e_value.tmp
# Paste
cd ${WORK_DIR}
paste gene_calls.tmp source.tmp accession.tmp function.tmp e_value.tmp > all.tmp
# header
echo -e "gene_callers_id\tsource\taccession\tfunction\te_value" > header.tmp
cat header.tmp all.tmp > ${WORK_DIR}/${MOTIF_NAME}_scan/gene-annot-CXXCH.tsv
rm *tmp
# Immport to contigs.db

source activate anvio-master; source ~/virtual-envs/anvio-master/bin/activate

anvi-import-functions -c $CONTIGS_DB -i ${WORK_DIR}/${MOTIF_NAME}_scan/gene-annot-CXXCH.tsv

deactivate; source deactivate
