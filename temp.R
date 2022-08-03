ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) 
    ncores
    if (is.na(ncores)) { 
        ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) 
    } 
    ncores
    if (is.na(ncores)) { 
    ncores<-2
    } 
    ncores