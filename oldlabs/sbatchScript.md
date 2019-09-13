# Sbatch script  
 
 ## How to start a sbatch script
 to start a sbatch script on uppmax you  type
 
    sbatch fileName.sbatch   
    
 
 
## An example of a sbatch file 
 
 
    #! /bin/bash -l
    #SBATCH -A b2015155
    #SBATCH -M milou
    #SBATCH -p core -n 16
    #SBATCH -t 1-00:00:00
    #SBATCH -J mRNAmappingPipeline_TRAP
    #SBATCH -e /proj/b2015155_nobackup/private/TRAP/reports/mRNAmappingPipeline_TRAP_SLURM_Job_id_%j.stderr.txt
    #SBATCH -o /proj/b2015155_nobackup/private/TRAP/reports/mRNAmappingPipeline_TRAP_SLURM_Job_id_%j.stdout.txt
    #SBATCH --mail-type=FAIL
    #SBATCH --mail-user=john.doe@scilifelab.se
    
    
    cd /pica/v10/b2015155_nobackup/private/TRAP
    # ==============================================================================
    # LOADING MODULES
    module load bioinfo-tools
    # modules for mapping
    module load star/2.4.1c
    module load samtools
    module load subread/1.5.0
    module use /proj/b2013006/sw/modules
    module load snakemake
    # modules for qc
    module load FastQC/0.11.2
    module load multiQC/v.0.2
    module load rseqc/2.6.2
    module load ucsc-utilities


    # Go to the working directory
    cd /pica/v10/b2015155_nobackup/private/TRAP
    
    
    # RUNNING SNAKEMAKE MAPPING SCRIPT WITH 16 cores
    snakemake -j 16 -s scripts/RNAmappingPipeline/SnakeMakeFiles/mappingReads.sm 

