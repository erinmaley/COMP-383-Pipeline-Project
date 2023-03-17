# COMP-383-Pipeline-Project
#
bowtie2, fastq-dump and SPAdes are needed for this code to run. Additionally, Biopython library is required.
# How to Download:
# bowtie2: 

$ wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-linux-x86_64.zip/download
or follow the instructions at http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

# fastq-dump:

follow the download and install instructions on https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft

# SPAdes:
follow the instructions on the download page, http://cab.spbu.ru/software/spades/#download.

# Biopython:

go into terminal and type pip!3 install biopython 

# Let's begin
In order to analyze the Genome Assembly of the human herpesvirus 5, Human cytomegalovirus(HCMV) in comparison to four Donors; Donor 1 (2dpi), Donor 1 (6dpi), Donor 3 (2 dpi), and Donor 3 (6dpi) we must first input the data from NCBI to our terminal

#
# 1
#

the links to the Donors are as follows: 

  Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360
  
  Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363
  
  Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374
  
  Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375
  
  
  First, follow the links provided respectfully and navigate by locating the SRR link under the "Runs:" catagory under the "Run" column, click on this link. Following this, locate the folder "Data access" between the "Reads" and "FAST/FASTQ" folder. Click on Data access. Scroll down and locate the highlighted blue link underneath the "SRA Normalized" table. Copy and keep this link for later.
  In terminal, use the following code to create a python shell:
  vim <INSERT_TITLE>.py
  in order to type press "I" 
  insert the following code in order to import the data we saved through the links earlier:
```python
  Donor 1 (2dpi)
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030')
  Donor 1 (6dpi)
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033')
  Donor 3 (2dpi)
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044')
  Donor 3 (6dpi)
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045')
```

press "esc" and then ":wq" to return to terminal, then:
python <SAME_TITLE_AS_EARLIER>.py 
to run the code

#
# 2
#

In order to compare strains we must formate the donor data, using the following command:

```python

os.system('fastq-dump -I --split-files SRR5660030')
os.system('fastq-dump -I --split-files SRR5660033')
os.system('fastq-dump -I --split-files SRR5660044')
os.system('fastq-dump -I --split-files SRR5660045')

```
next, using Bowtie2, we must create and index for HCMV, in the following command:
```python
import os
import Bio
from Bio import Entrez

Entrez.email = "email@mail.com"

handle = Entrez.efetch(db="nucleotide",id='NC_006273.2',rettype='fasta')

fasta = handle.read()

handle.close()

with open ('HCMV.fasta','w') as f:

  f.write(fasta)
```
this will read in and format HCMV into a fasta file
next we perform bowtie2 in order to compare the overlaps to the donors genome and HCMV:

```python
os.system('bowtie2-build HCMV.fasta HCMV)
#this creates an index in order to perform bowtie 2

os.system('bowtie2 --quiet -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S HCMV30mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S HCMV33mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660044_1.fastq -2 SRR5660044_2.fastq -S HCMV44mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -S HCMV45mapped.sam')
```
this is used to read in the following donor data and compare it too HCMV


```python
#Getting initial read counts for each sample
initial_read_counts = {}
samples_list = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
for s in samples_list:
   with open(f'{s}_1.fastq') as file1:
       initial_read_counts[s] = sum([1 for line in file1]) / 4
   print(f'{s} has {initial_read_counts[s]:,} read pairs before filtering')
```
```python
# Filtering out non-mapped reads and writing mapped reads to new fastq files
for s in samples_list:
   output_sam = f'{s}_mapped.sam'
   os.system(f'bowtie2 --quiet -x HCMV -1 {s}_1.fastq -2 {s}_2.fastq -S {output_sam}')

   post_filter_counts = 0
   with open(output_sam) as samfile, open(f'{s}_HCMV.fastq', 'w') as fastqfile:
       for line in samfile:
           if line.startswith('@'):
               continue
           parts = line.split('\t')
           flag_val = int(parts[1])

           if flag_val & 4 == 0:
               post_filter_counts += 1
               qname_val = parts[0]
               seq_val = parts[9]
               qual_val = parts[10]
               fastqfile.write(f"@{qname_val}\n{seq_val}\n+\n{qual_val}\n")
              
```
```python
   # Writing filtered read counts to a log file
   with open('log.txt', 'a') as logfile:
       logfile.write(f'{s} has {initial_counts[s]:,} read pairs before filtering and {post_filter_counts:,} read pairs after filtering.\n')

```

#
# 3
#

```python
# import the os module
import os

# create a dictionary containing sample names and their corresponding fastq file names
samples = {
    'SRR5660030': ('SRR5660030_HCMV.fastq'),
    'SRR5660033': ('SRR5660033_HCMV.fastq'),
    'SRR5660044': ('SRR5660044_HCMV.fastq'),
    'SRR5660045': ('SRR5660045_HCMV.fastq'),
}

# initialize an empty string to store the SPAdes input command
spades_input = ''

# iterate through the samples dictionary and construct the SPAdes input command
for index, (<sample_name>, fq) in enumerate(samples.items(), start=1):
    spades_input += f'--s{index} {fq} '

# construct the SPAdes command with the input command and other parameters
spades_command = f'spades.py -k 77,99,127 -t 4 --only-assembler {spades_input.strip()} -o HCMV_SRR_assembly'

# execute the SPAdes command using the os.system() method
os.system(spades_command)

# append the SPAdes command to a log file
with open('log.txt', 'a') as log_file:
    log_file.write(f'SPAdes command: {spades_command}\n')
```
This step should take awhile
#
# 4
#
```python
# import the SeqIO module from the Bio library
from Bio import SeqIO

# set the output file name from the assembly
output_file = 'contigs.fasta'

# initialize variables to keep track of the number of contigs greater than 1000 bp and the total assembly length
contigs_gt_1000 = 0
assembly_length = 0

# iterate through the contigs in the output file from the assembly and calculate the length of each contig
for record in SeqIO.parse(output_file, 'fasta'):
    contig_length = len(record.seq)

    # if the length of a contig is greater than 1000 bp, increment the count of contigs greater than 1000 bp and add the length of the contig to the assembly length
    if contig_length > 1000:
        contigs_gt_1000 += 1
        assembly_length += contig_length

# append the number of contigs greater than 1000 bp and the total assembly length to a log file
with open('log.txt', 'a') as log_file:
    log_file.write(f'There are {contigs_gt_1000} contigs > 1000 bp in the assembly.\n')
    log_file.write(f'There are {assembly_length} bp in the assembly.\n')
```

#
# 5
#
```python
from Bio.Blast import NCBIXML

# Open the BLAST result file in read mode and parse it into a list of Blast records
with open("blast_results.xml", "r") as result_handle:
    blast_records = list(NCBIXML.parse(result_handle))

# Open a new file in write mode to store the top 10 matches and write the header
with open("<file_name_top>.txt", "w") as log_file:
    log_file.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")

    # Loop through each Blast record
    for record in blast_records:
        
        # Loop through the top 10 alignments for each record
        for alignment in record.alignments[:10]:  
            # Extract the first HSP for each alignment
            hsp = alignment.hsps[0]  

            # Format the output line with the required fields and write it to the log file
            output_line = f"{alignment.accession}\t{hsp.identities * 100 / hsp.align_length:.2f}\t{hsp.align_length}\t{hsp.query_start}\t{hsp.query_end}\t{hsp.sbjct_start}\t{hsp.sbjct_end}\t{hsp.bits}\t{hsp.expect}\t{alignment.title}\n"
            log_file.write(output_line)
```
# Done
