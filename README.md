# Usadel for nanowires

This code solves the Usadel equation in a superconductor-nanowire-superconductor (SNS) junction. An in plane exhange field and spin-orbit terms are also impelemented on the code. There are two run scripts, one for the supercurrent through the junction and one for the density of states (DOS). For both, the supercurrent and the DOS, there is also a corresponding plot script.

The code is operated using run_script_supercurrent.sh for the supercurrent and run_script_dos.sh for the DOS and all the plotting is done by using the corresponding plot scripts. Explanations for available parameters are found in the scripts.

Some results obtained with the code are presented in INSERT REFERENCE. See [On citations](##on\ citations).

## Installation guide for the dependencies

### Setuptools

If you do not have setuptools installed install it by running

```
sudo apt-get install python-setuptools
```

### fimport

Download the fimport from https://pypi.python.org/pypi/fimport and extract the tar.gz file to a folder for example in /home/Programs/. The choice of folder doesn't matter. Move to the folder in terminal with

```
cd /home/Programs/fimport-version_number/
```

and install fimport with

```
sudo python setup.py install
```

### gfortran

If you do not have a fortran compiler, it is recommended to install gfortran. This can be done with

```
sudo apt-get install gfortran-4.8
```

or by searching "gfortran" in Ubuntu software center and installing it from there.

### scikits.bvp_solver

The easiest way to install scikits.bvp_solver is by downloading the python egg file from https://pypi.python.org/pypi/scikits.bvp_solver, moving into the directory where it is, and running 

```
sudo easy_install scikits.bvp_solver-1.1-py2.7-linux-i686.egg
```

from the terminal. You can also download the tar.gz file and compile the scikits.bvp_solver in the same way as you did with fimport.

If you installed scikits.bvp_solver by using the python egg, you might get a warning

/usr/lib/python2.7/dist-packages/pkg_resources.py:1031: UserWarning: /home/username/.python-eggs is writable by group/others and vulnerable to attack when used with get_resource_filename. Consider a more secure location (set with .set_extraction_path or the PYTHON_EGG_CACHE environment variable).
  warnings.warn(msg, UserWarning)

when running the code. You can get rid of the warning by changing the ~/.python-eggs directory not writeable by everyone with

```
chmod g-wx,o-wx ~/.python-eggs
```

in you home folder.

Now you should have all the dependencies to run the code.

## Example

You can test the code after installation by running the run_script_dos_test.sh from the terminal while being in the main directory of the code

```
bash Test/run_script_dos_test.sh
```

The script will calculate the density of states for exchange fields h = 0.00, 8.00, 16.00, and 24.00 with the SO term A = (0.00, 0.00, 2.00, 0.00). The values of the rest of the parameters can be found in the scripts. You can then plot the obtained curves with the plot_script_dos_test.sh by running

```
bash Test/plot_script_dos_test.sh
```

in the terminal. The figure is by default stored in a Test/DOS_figures subfolder.

You can run similar test for the supercurrent with 

```
bash Test/run_script_supercurrent_test.sh
```

and then

```
bash Test/plot_script_supercurrent_test.sh
```
in the main directory of the code.

The first script will calculate the supercurrent for the SO term A = (0, 0, alpha sigma_x, 0) with alpha=0.00, 1.00, 2.00 for the exchange fields from 0.00 to 80.00. The second script will plot the results and store the figure in the Test/Figures subfolder.

Note that the actual solving of the equations takes quite a long time even though the matsubara_sum parameter is set to 50 in the run script, meaning that the results aren't yet that accurate. For accurate results use the default value matsubara_sum=0 which adjusts the accuracy depending on the chosen temperature. Also the tolerance parameter is set to 10e-4 for faster results. I would recommend using tolerance=10e-6 in actual calculations.

## On citations

If you use this code in academic publications, I encourage you to cite it appropriately, for example in BibTeX format

```
@Article{
    author = {J. Arjorant and T.T. Heikkilä}
    title = {Intrinsic spin-orbit interaction in diffusive normal wire Josephson weak links: supercurrent and density of states},
    journal = {Phys. Rev. B},
    volume = {},
    pages = {},
    year = {2016},
    note = {Code available at https://github.com/wompo/Usadel-for-nanowires/}
}
```

## Note

The example run scripts also clean up the solution files after they have finished by removing the created solutions subdirectories. The actual run scripts in the main folder of the code to not remove the solutions automatically since they can be used later as initial guesses for future calculations. Please do note that the solution files tend to build up and take up a lot of disk space. Therefore manual cleaning up is sometimes necessary. All the solutions are stored in subfolders defined in the run scripts.

The code has been tested with

- python 2.7.6
- numpy 1.8.2
- gfortran 4.8
- fimport 0.2
- scikits.bvp_solver 1.1
- setuptools 3.3
