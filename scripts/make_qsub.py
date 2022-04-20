import glob

for file in glob.glob("*_1.fastq.gz"):
    filename = file.split("_1.fastq.gz")
    job_name = filename[0] + ".pbs"
    f = open(job_name, "w")

    f.write("#!/bin/bash\n\n#PBS -l select=1:ncpus=16:mem=16gb\n#PBS -l walltime=12:00:00\n")
    f.write("\nmodule load bowtie/2.3.5.1\nmodule load metaphlan/3.0.10\n")
    f.write("cd /srv/scratch/vafaeelab/AbhishekVijayan/sample_size_estimation/PRJEB43119/data\n\n")

    # metaphlan ERR6230245_PRIMM0513_1.fastq.gz,ERR6230245_PRIMM0513_2.fastq.gz --input_type fastq --bowtie2out ERR6230245_PRIMM0513.bowtie2.bz2 --nproc 4 > ERR6230245_PRIMM0513_profile.txt

    metaphlan_command = ("metaphlan " + file + "," + filename[0] + "_2.fastq.gz" +
        " --input_type fastq --nproc 16 --bowtie2out " +  filename[0] + ".bowtie2.bz2 > " 
        + filename[0] + "_profile.txt")

    f.write(metaphlan_command)

    f.close()