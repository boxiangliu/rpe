Author: Mike Gloudemans

Existing RPE GWAS and eQTL files were formatted prepared for colocalization analysis with the following command:

`bash rpe_prep.sh`

I then ran the FINEMAP colocalization pipeline on these files. 

In this folder, I have placed the colocalization pipeline config file `settings_used.config` along with the `git log` at
the time this pipeline was run (`git_status.txt`). If it is ever necessary to rerun the pipeline in exact the same way, we could do so by
checking out the appropriate version from the git repository, and running the old pipeline using the included config file.
