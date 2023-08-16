
Aether uses containers to create storage space for its elements. In 
this scenario, containers are used in the means to manage and organize 
output data created during simulation. Due to its efficient structure, 
containers allow for easy storage and accessibility which makes larger 
information easier to handle. This is useful for the visualization and 
analysis of the data.

In the OutputContainer class, the vector AllOutputContainers is initialized
to store all the OutputContainers for different output types. The container
is filled with variables and data specified by the output type.

An important trait of containers is their flexible data storage as it adjusts 
to the size of data given. Their dynamic sizing abilities creates a unique 
container for each type of output. 

As an example: 
In the OutputContainer class, if the output_type is ions or states, its 
respective container will store the "Bulk Ion" data (via store_variable).
    if (type_output == "ions" || type_output == "states")
        AllOutputContainers[iOutput].store_variable("Bulk Ion " +
                                                    ions.temperature_name,
                                                    ions.temperature_unit,
                                                    ions.temperature_scgc);

Data is outputted from the containers using "write" method. They can write data 
to two types of files, netCDF or binary file. For example: output_binary_3d 
function. Data can also be read from files and stored in containers using 
the "read" method. They can read from netCDF and binary file formats. For 
example: input_binary_3d. By using these file-specific methods, Aether
employs a more efficient way of handling data.
