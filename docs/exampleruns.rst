=========================================
A guide to running the model
=========================================
Here is a quick guide to running the model.

Paths
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The example scripts define a few paths, required to run the model:
- src_dir: This is the directory to the source code. Please change accordingly
- save_path: This is where data from the model runs will be saved.
  I suggest making separate directories for the init runs and the SS runs (and for each WC value).
- fig_path: This is where automatically generated figures will be saved.


Saving data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
There are two ways to save data.
The first is by setting save_bool = true.
This will save snapshots of the streamfunction
at periods specified via the save_every variable.


The second option for saving data is by setting save_last=true.
This will save a complete snapshot of the upper and lower streamfunction at the
end of the model run, so that the snapshot can be used for initial conditions
for future runs.


Plotting
~~~~~~~~~~~~~~~~~~~~~~~~~~~
A plotting function is currently defined in the output_fcns.jl file.
How often plots are made is set by the plot_every parameter.

Running the model
~~~~~~~~~~~~~~~~~~~~~~~~~
Run the BT_test.jl file to run the model.
All adjustable parameters can be changed in that file.


