# ABC - SLRUM - TWITTER-SPREAD

This Folder `abcdir` is a working folder yet not well organized but that aims at storing all scipt needed to run and analyse the model available in `twitter-spread` usings ABC. 

I will keep some note here on the name of the folders: usually VF or TF mean "True False/ Vrai Faux". This generally note that the abc has been run with the score computed against both distribution of false and true tweet. `testneutralRumors` means `Rumors` _ie_ we test 3Dcascade on the distribution of RUMORS and not on all cascades.


## Content:

All fine `abc*.R` (like: `abcRandom.R`, `abcRumors.R`,...) are script able to run ABC using MPI on their own. They generate priors, run the experiment using `parLapply` function from `Rmpi` package. They create output folder where the scrore and the parameters of the simulations will be stored.

This is ok for relatively quick model that all last around the same time. This is not the case with the 3D model that can sometime generate really really slow simulations, which will then ask SLURM for hours when just one node will need this hours. To try to optimize this we wrote a "generateParameters.R"  

- file in the folder `slurm-mpi/` are script used to run bunch of experiment via mpi and SLURM.

## file `generateParameters.R`
The idea behind this script is to generate a huge set of parameter that can be split given their expected running time. The times used in the current version of this folder have been calculated by running bunch of test and using a log-linear fitting. Equation and  should be somewhere in the scripts here .


To generate the folder with all parametesr use something like:

```bash
Rscript generateParameters.R numfolder_stored_parameters $folder_stored_parameters
```
This will be read using  the script `readParameters.R` and run using `Cascade3D`.  `readParameters.R` is wrapped within `slurm-mpi/Rmpiscript-mn4_splitted_large` and `slurm-mpi/Rmpiscript-mn4_splitted_small` scripts that wrap the call within SLURM scheduling system and MPI.

Then you can use the `slurm-mpi/send_batch_simu` file that look like:

```bash
folder_stored_parameterss=nonneutralRumors
outputSimulation=nonneutralRumorsVF
for i in {3..5};
do
    echo "running $folder_stored_parameterss/$folder_stored_parameterss$i"
    sbatch Rmpiscript-mn4_splitted_large $folder_stored_parameterss/$folder_stored_parameterss$i $outputSimulation "large"
    sbatch Rmpiscript-mn4_splitted_small $folder_stored_parameterss/$folder_stored_parameterss$i $outputSimulation "small"
done
```

