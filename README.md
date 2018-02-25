# cc-abide

# repo structure overview

`merge.R` merges data from `data/abide/` and `data/pcp/`, and outputs to `data/`.

# `data/pcp/`
Some data in `data/pcp/` are loaded and saved as large number of separate files
(one for each of 1035 patients). To speed up and simplify file synchronization
they are packed as `data/pcp/freesurfer_raw_data.tar.gz`. If there's a need to
read these individual files you can easily unpack the archive, `.gitignore` ignores 
the unpacked contents.

Downloading PCP data: `data/pcp/download.sh` takes `data/pcp/FILE_ID.txt` as input
(generated by `merge.R` from `data/pcp/summary.csv`). You have to manually
edit `download.sh` and set the proper URL template, depending what you want,
see [PCP download instructions](http://preprocessed-connectomes-project.org/abide/download.html).

# `merge.R`
I tried to keep the original varaible names from ABIDE and PCP as much
as possible, I left a `# VARIABLE RENAME` comment line wherever those
names are changed for simplification or consistency.

All merged tables generated by `merge.R` are saved as .csv in `data/` and may
be loaded to R directly, without the need to run the whole script.
`A.csv` is the last table and includes data from all others.
