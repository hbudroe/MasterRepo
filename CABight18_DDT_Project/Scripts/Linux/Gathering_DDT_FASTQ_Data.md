## Small file for initially obtaining the raw ddt fastq files from mars8180 teaching cluster & scripts/process through qiime2

### Transfer data
Login to xfer node
Transfer zipped file from my ddt-project/data folder in the teaching cluster (originally copied from instructor data)
``` bash
scp -r -p /home/hmb25721/ddt-project/data/ddt-raw-fastq.tar.gz sapelo2.gacrc.uga.edu:
```
### Check on data
Login to sapelo2 
Move (not copy) file from home to the ddt_project/data folder in my home dir
Unzip the tar file
Check a random sample to make sure it worked
``` bash
mv ddt-raw-fastq.tar.gz ddt_project/data/
tar -xvzf ddt-raw-fastq.tar.gz
zcat DDT.14.2_S1358_L001_R2_001.fastq.gz | head -n 8 #or any sample name is fine #should print 4 lines (sequencing info, sequence, +, phred score
```
## Run scripts:
### 01- "import" to qiime -- aka convert the original fastq file into a qiime artifact .qza
  Can see what that looks like 
```
cat ddt_project/results/01-convert-to-qiime.qza  | head -n 8
```
### 02- qiime summarize -- summarize our data quality into a visual .qzv
  pain in the butt, but login to xfer node to transfer the .qzv to local desktop so I can
  drag/drop file into https://view.qiime2.org/
  BUT 5/28/25 I spent 2 hours trying to scp and kept getting "No such file or directory" even though I can ls -lsa to the exact path & have write permissions :/ Gave up temporarily and used globus to download, but     come back here!!!
### 03 - trim adaptors with cutadapt & summarize post-primer trimming quality as a .qzv
### 04 - run dada2 clustering/denoising 
  trunc-len-f/r to decide read truncating based on quality visualization from 03 summary
  mars 8180 alejandro did 290 / 290 ; I think maybe should reduce -r to ~ 215ish since 290 seems too far into the read with not great quality, leave -f as 290 for now
### 05 - assign taxonomy using blast in qiime2
  Use new modified silva database from Tiago Mar-10-2025 currently stored in /project folder (can only be accessed via xfer node [txfer is teaching transfer node silly)
  Database paths: /project/hmblab/lab_shared/silva-taxonomy/ref-dbs/SILVA138-nr99_fixed_strings_custom_seq_sequences-Mar-10-2025.qza &  SILVA138-nr99_fixed_strings_custom_seq_taxonomy-Mar-10-2025.qza
  Also 05.1 - assign taxonomy with 80% similarity (blast-80-1)
  ```
#while logged into xfer node:
cp /project/hmblab/lab_shared/silva-taxonomy/ref-dbs/SILVA138-nr99_fixed_strings_custom_seq_sequences-Mar-10-2025.qza /home/hmb25721/ddt_project/databases
cp /project/hmblab/lab_shared/silva-taxonomy/ref-dbs/SILVA138-nr99_fixed_strings_custom_seq_taxonomy-Mar-10-2025.qza /home/hmb25721/ddt_project/databases
```
  Then run script 05 using the new database paths as input
Below are some comments that included on my original script that I think are the reason for the immediate job failure but hepful for my brain to undersatnd
  ```
qiime feature-classifier classify-consensus-blast \
  --i-query ${INPUT} \   #input sequence table
  --i-reference-taxonomy ${REFTAX} \ #reference tax labels
  --i-reference-reads ${REFSEQ} \ #reference sequences
  --p-maxaccepts 1 \ #max number hits to keep (default 10)
  --p-perc-identity 0.90 \ #reject of %identity lower
  --o-classification ${CLASS} \  #output taxonomy classifications
  --o-search-results ${SEARCH} \ #top hits for each query
  --p-num-threads 12
# --i = inputs ; --p = parameters ; --o = outputs
```
# 06 - build phylogenetic tree
  For this one, there was no reference script for mars8180 (we just downloaded an existing tree from instructor data), so try bulding own & troubleshoot as needed
  Troubleshooting : failed the 1st time --> --p-num-threads change to --p-n-threads
  But it worked the second time! woohoo
```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${INPUT} \
  --o-alignment ${ALIGNMENT} \
  --o-masked-alignment ${MASKEDALIGNMENT} \
  --o-tree ${UNROOTEDTREE} \
  --o-rooted-tree ${ROOTEDTREE} \
  --output-dir ${OUTPUT} \
  --p-n-threads 12
```
In the interim, download the asv/tax/metadata
```
#From the teaching cluster download the metadata
scp hmb25721@txfer.gacrc.uga.edu:/work/mars8180/instructor_data/metabarcoding-datasets/ddt-project/metadata/2025-01-03-ddt-metadata.csv .
#From my home directory download the asv feature table & assigned taxonomy
scp hmb25721@xfer.gacrc.uga.edu:/home/hmb25721/ddt_project/results/04-dada2-feature-table.qza .
scp hmb25721@xfer.gacrc.uga.edu:/home/hmb25721/ddt_project/results/05-taxonomy-blast-90-1.qza .
```

# 06.1 - build phylogenetic tree using qiime2::iqtree
  There are some (possible) issues with the phylogeny built in the fasttree tree... For ex, Gammanema doesn't appear placed with the other genus in its family (Selachinematidae:Halichoanolaimus) -- this could be an issue with taxonomy assignment or the tree. Rebuilding the phylotree using a better algorithm to assess & for future phylogenetic analyses.
  iqtree requires a multiple sequence alignment as input. Thankfully, the fasttree command I ran previously already output an alignment (alongside the tree) using mafft so I can use that as input instead of running de novo msa with mafft.
```
qiime phylogeny iqtree \
        --i-alignment ${INPUT}\
        --p-lbp 10000 \
        --o-tree ${OUTPUT}
```
