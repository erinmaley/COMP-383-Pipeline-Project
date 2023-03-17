# COMP-383-Pipeline-Project
#In order to analyze the Genome Assembly of the human herpesvirus 5, Human cytomegalovirus(HCMV) in comparison to four Donors; Donor 1 (2dpi), Donor 1 (6dpi), Donor 3 (2 dpi), and Donor 3 (6dpi) we must first input the data from NCBI to our terminal

-1-


#the links to the Donors are as follows: 

  #Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360
  
  #Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363
  
  #Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374
  
  #Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375
  
  
  #First, follow the links provided respectfully and navigate by locating the SRR link under the "Runs:" catagory under the "Run" column, click on this link. Following this, locate the folder "Data access" between the "Reads" and "FAST/FASTQ" folder. Click on Data access. Scroll down and locate the highlighted blue link underneath the "SRA Normalized" table. Copy and keep this link for later.
  #In terminal, use the following code to create a python shell:
  vim <INSERT_TITLE>.py
  #in order to type press "I" 
  #insert the following code in order to import the data we saved through the links earlier:
  #Donor 1 (2dpi)
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030')
  #Donor 1 (6dpi)
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033')
  #Donor 3 (2dpi)
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044')
  #Donor 3 (6dpi)
os.system('wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045')

#press "esc" and then ":wq" to return to terminal, then:
python <SAME_TITLE_AS_EARLIER>.py 
#to run the code

-2-

#In order to compare strains we must formate the donor data, using the following command:

#os.system('fastq-dump -I --split-files SRR5660030')
#os.system('fastq-dump -I --split-files SRR5660033')
#os.system('fastq-dump -I --split-files SRR5660044')
#os.system('fastq-dump -I --split-files SRR5660045')

#next, using Bowtie2, we must create and index for HCMV, in the following command:

import os
import Bio
from Bio import Entrez

Entrez.email = "email@mail.com"

handle = Entrez.efetch(db="nucleotide",id='NC_006273.2',rettype='fasta')

fasta = handle.read()

handle.close()

with open ('HCMV.fasta','w') as f:

  f.write(fasta)
  
#this will read in and format HCMV into a fasta file
#next we perform bowtie2 in order to compare the overlaps to the donors genome and HCMV:

os.system('bowtie2-build HCMV.fasta HCMV)
#this creates an index in order to perform bowtie 2

os.system('bowtie2 --quiet -x HCMV -1 SRR5660030_1.fastq -2 SRR5660030_2.fastq -S HCMV30mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660033_1.fastq -2 SRR5660033_2.fastq -S HCMV33mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660044_1.fastq -2 SRR5660044_2.fastq -S HCMV44mapped.sam')
os.system('bowtie2 --quiet -x HCMV -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -S HCMV45mapped.sam')

#this is used to read in the following donor data and compare it too HCMV

