These files are provided 'as is' to provide a record of the figures, data and analysis presented in the associated paper. They are messy. Scripts are not designed to work on your machine - although you should be able to run the figure construction scripts in R if you have the necessary packages installed. 

Any questions, please contact:
liam.shaw@biology.ox.ac.uk
liam.philip.shaw@gmail.com

Author: Liam Shaw, 05/09/2022
Updated: 28/03/2023

Figshare does not have a folder structure, so I use prefixes to organise the files. 

## PREPROCESSING 
`preprocessing_trim-and-assemble.sh` - cluster array script used to download reads from NCBI then de novo assemble
`preprocessing_map-with-snippy.sh` - cluster array script used to map reads (after NCBI download) to a reference genome
`preprocessing_snippy-core.sh` - cluster script (single command) used to get variants from snippy mapping, with masked positions

## ASSEMBLY DATA
`assemblies.zip` - de novo assemblies for all isolates analysed

## OTHER DATA FILES 
`data_mask.bed` - BED file of positions masked (for snippy-core)
`data_metadata.csv` - metadata for isolates
`data_shen-mlst.tsv` - ST assignments for isolates, from de novo assemblies and mlst (Torsten Seeman, https://github.com/tseemann/mlst)
`data_snippy-core-results*` - output files from snippy-core
`data_KU761327.1.gff3.tsv` - annotations for the reference IncX4 plasmid KU761327.1


## FIGURE CONSTRUCTION (R scripts)
`script_Figure-5A.R` -- makes Sankey diagram, panel A 
`script_Figure-5B.R` - panel B, minor variant results and gene diagram
`script_Figure-5C.R` - IncX4 phylogeny, panel C
`script_Figure5.R` - makes combined plot (Figure 5 of manuscript). 
Note that these figures were subsequently edited in Inkscape for publication. 

## SUPPLEMENTARY FIGURE CONSTRUCTION (R SCRIPTS)
`script_supplementary-IncX4-phylogeny.R` - makes emergence phylogeny showing PV2 clade within PV3. Rooted to reference.
`script_supplementary-IncX4-PV3-time.R` - makes bar charts showing IncX4 and PV3 association over time. 


