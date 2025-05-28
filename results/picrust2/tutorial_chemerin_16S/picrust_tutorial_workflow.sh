##### tutorial for picrust2

#### NOTE you must have activated the picrust2 virtualenv
conda activate picrust2

#### Getting started

### get data (done)

# wget http://kronos.pharmacology.dal.ca/public_files/picrust/picrust2_tutorial_files/chemerin_16S.zip
# unzip chemerin_16S.zip
cd tutorial_chemerin_16S

### Looking at input files required

ls -1
# metadata.tsv
# seqs.fna
# table.biom

less seqs.fna
# >86681c8e9c64b6683071cf185f9e3419
# CGGGCTCAAACGGGGAGTGACTCAGGCAGAGACGCCTGTTTCCCACGGGACACTCCCCGAGGTGCTGCATGGTTGTCGTC
# AGCTCGTGCCGTGAGGTGTCGGCTTAAGTGCCATAACGAGCGCAACCCCCGCGTGCAGTTGCTAACAGATAACGCTGAGG
# ACTCTGCACGGACTGCCGGCGCAAGCCGCGAGGAAGGCGGGGATGACGTCAAATCAGCACGGCCCTTACGTCCGGGGCGA
# CACACGTGTTACAATGGGTGGTACAGCGGGAAGCCAGGCGGCGACGCCGAGCGGAACCCGAAATCCACTCTCAGTTCGGA
# TCGGAGTCTGCAACCCGACTCCGTGAAGCTGGATTCGCTAGTAATCGCGCATCAGCCATGGCGCGGTGAATACGTTCCCG
# >b7b02ffee9c9862f822c7fa87f055090
# AGGTCTTGACATCCAGTGCAAACCTAAGAGATTAGGTGTTCCCTTCGGGGACGCTGAGACAGGTGGTGCATGGCTGTCGT
# ...

    # file must not contain raw (non denoised/unclustered) reads
    # header lines (with >) must only contain sequence ID
    # seqs must be positive strand of 16s rRNA

biom head -i table.biom
    
    # first col- IDs of FASTA 
    # other cols- sample IDs
    # input table must be ABSOLUTE ABUNDANCE
        # NOT RELATIVE ABUNDANCE
        # can be rarefied table

biom summarize-table -i table.biom

    # use to check + remove any low-depth samples

##### run the pipeline

#### run whole pipeline
picrust2_pipeline.py -s study_seqs.fna -i study_seqs.biom -o picrust2_out_pipeline_split -p 11

   # 11 cores used 

#### run pipeline step by step

### mkdir for output files

mkdir picrust2_out_pipeline
cd picrust2_out_pipeline

### Place reads into reference tree

# done by place_seqs.py
    # aligns seqs 
        # using HMMER
    # finds most likely placements in reference tree
        # using EPA-NG or SEPP
    # outputs treefile w ASVs as tips
        # using GAPPA

place_seqs.py -s ../seqs.fna -o out.tre -p 11 \
              --intermediate intermediate/place_seqs
    # the 11 in the following commands is num cores to use
    # last flag is where intermediate files go (idc)

### Hidden-state prediction of gene families, hsp

hsp.py -i 16S -t out.tre -o marker_predicted_and_nsti.tsv.gz -p 11 -n

hsp.py -i EC -t out.tre -o EC_predicted.tsv.gz -p 11

# predicts missing genome for each ASV
    # preds copy num of gene families for each ASV

# outputs
zless -S marker_predicted_and_nsti.tsv.gz
    # NSTI- how closely related taxa are to reference genomes
        # average phylo distance (substitutions per site) to nearest relative
    # less but for gzipped files
        # -S flag turns off line wrapping **

zless -S EC_predicted.tsv.gz
    # predicted copy num of all enzyme classification nums 
        # for each ASV
    # like dewey decimal for enzymes- 1.1.1.1 = alcohol dehydrogenase
    # what's used to predict metacyc pathway levels

### Generate metagenome predictions

metagenome_pipeline.py -i ../table.biom -m marker_predicted_and_nsti.tsv.gz -f EC_predicted.tsv.gz \
                       -o EC_metagenome_out --strat_out
    # theory
        # control for variation in 16S copy numbers across orgs
            # divide read depth per ASV by predicted 16S copy numbers
        # multiply normalized^ ASV read depth per sample by predicted gene family copy numbers per ASV
    # code
        # --strat_out -> stratified output file
            # = long format, shows how ASVs contribute ECs

# output
zless -S EC_metagenome_out/pred_metagenome_unstrat.tsv.gz
    # ECs in each sample

zless -S EC_metagenome_out/pred_metagenome_contrib.tsv.gz
    # taxon_rel_abun
        # sum of all taxa abundances per sample is 100).

    # genome_function_count
        # Predicted copy num of this fn per taxon.
        # think of fn like enzyme gene

    # taxon_function_abun
        # Multiply "taxon_abun" by "genome_function_count".

    # taxon_rel_function_abun 
        # Multiply "taxon_rel_abun" by "genome_function_count".

    # norm_taxon_function_contrib
        # Proportional rel abun of the taxon_function_abun per function and sample
        # I.e what proportion of the specified function is contributed by that particular taxon in the sample.

# to convert to legacy format for burrito/mimosa
# convert_table.py EC_metagenome_out/pred_metagenome_contrib.tsv.gz \
                #  -c contrib_to_legacy \
                #  -o EC_metagenome_out/pred_metagenome_contrib.legacy.tsv.gz

### Pathway-level inference
pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_contrib.tsv.gz \
                    -o pathways_out -p 11

# what it does
    # regroups EC nums to MetaCyc reactions ie finds what the enzymes do
    # infers Metacyc pathways based on these
    # returns abundance of pathways identified

# output
    # how much each ASV is contributing to community-wide pathway abundance
        # could also do --per_sequence_contrib option

### Add functional descriptions
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz

# gives a description of each fnal ID in output tables
