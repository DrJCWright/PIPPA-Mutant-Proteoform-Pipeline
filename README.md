**PIPPA Mutant Proteoforms Pipeline**

_Dr James C Wright
2025 The Institute of Cancer Research_

MIT License

A Pipeline for Mutant Proteoform Discovery, Quantification, and Annotation

This is a complex NextFlow pipeline using custom tools combined with OpenMS tools.
We recommend that you familiarise yourself with both NextFlow and OpenMS before running PIPPA:
    https://openms.de/
    https://www.nextflow.io/

**INSTALLATION:**

0. Install Nextflow on your system:  https://www.nextflow.io/

1. PIPPA uses the bioconda OpenMS release and OpenMS 3rd Party Tools for many functions. 
    You will need to run the pipeline on a system with conda available.
    We recommend installing Mamba for fast building of enviroments:
        https://github.com/mamba-org/mamba


2. Build OpenMS Conda Enviroment using the openms-conda-enviroment.yml from PIPPA root directory run:

    >mamba env create -f openms-conda-enviroment.yml -p ./env/openms

    If on a Slurm HPC system you can sbatch the install_openms.sh script to build the required enviroment.
    If not using mamba you can substitute "mamba" with "conda" or you preferred conda tool.

3. If you want to use MSFragger you will need to download it from here: 
    https://github.com/Nesvilab/MSFragger/wiki/Preparing-MSFragger#Downloading-MSFragger 

    Once downloaded you can specify the path to the MSFragger Jar file in the nextflow.config as the param "msfragger_path"

4. Check params for the MSGFplus exacutable (installed as part of the OpenMS conda enviroment) "msgf_path"

5. Download reference sequences from GENCODE (https://www.gencodegenes.org/) and place them in the "data" directory:
        a. gencode.*.pc_translations.fasta
        b. gencode.*.pc_transcripts.fasta
    Update the following params in nextflow.config:
        cosmic_aa_reference =
        cosmic_na_reference =
        reference_fasta     =

    Additionally you can use or append UniProt or other canonical protein fasta to the reference_fasta. (This is used to label canonical peptides and proteins vs non-canonical or variants )

6. Download COSMIC mutations tsv:
    Login or create and account for COSMIC here: https://cancer.sanger.ac.uk/cosmic/login
    Navigate to data downloads and download the relevant "Genome Screens Mutant" Tsv file (make sure genome version matches the GENCODE reference)
    Put the file in the data directory and update nexflow.config param "mutation_table"


**RUNNING PIPPA MUTANT PROTEOFORMS:**

1. Prep Sample Information
    create a tsv file with deatils of the TMT experiment you are processing. The TSV should have the following columns:
    ExperimentID, TMTLabel, SampleID

    (A basic example file is provided in the data directory)
    Check the param "sample_info" points to the table this can be done at run time or in submission script using --sample_info='./path_to/tsv_file.tsv'

2. Prep input spectra:
    Convert spectra to standard mzML format and put the mzML files in the mzml directory.

    NB. The default location is the "mzml" directory in the root of PIPPA however you can specify a custom location using the param "mzml_dir" (i.e. at runtime or in slurm submission script using --mzml_dir='./mzml/' )

    If spectra are fractioinated PIPPA can merge fractions into a single search using the '--mzml_merge_pattern' this is a regex pattern that will group spectra together. (i.e '.*' will treat each mzML file as a seperate sample, 'sample_name' will group all files matching sample_name, and for more advanced combinations regex patterns such as 'sample_name_run\d+' will merge mzMLs but only where the run number matches. )

3. Prep databases:
    PIPPA will build proteoform databases and add decoys using mutation information and gencode references, however, it is recommended you include additional reference fasta database(s) in the db directory. We suggest at minimum putting a copy of gencode.*.pc_translations.fasta in the db directory.

    NB. As above you can override default database directory using param --database_dir 

    If you are rerunning data and have already build the proteoform databases then you can skip the build proteoform database script using the param:
        --build_mutant_db='false'
    Make sure that if this is the case that the mutation info tsv generated during database building is set using the param:
        --mutation_info='./path/to/mutation_info.tsv'

4. Prep run script and check search params:
    An example run script can be found in 'run_scripts' directory it is designed for Slurm but can easily be adapted for local systems by dropping the 'slurm' profile. 

    Setting MS search params can be done by overiding them in the run script or creating a custom Profile in the nextflow.config.
    A standard TMT 10 plex MS3 profile is already provided.

5. you can now run PIPPA using the modified run_script or using the command:
    >nextflow run main.nf -profile slurm,TMT10_MS3

    Results will be written to the directory specified by '--out_dir' there will be multiple subdirectories but the main final output tables with PSMs, peptides, proteins, and proteoforms will be in 'final_results'

    NB. Avoid using the same results folder for multiple runs as previous files will get overwritten or reused. Either rename the results folder when PIPPA has completed or better point 'out_dir' to different results folder for each run.

    There are many parameters in the nextflow.config that can be configured. Most are described in the OpenMS documentation.

    




