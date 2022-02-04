// Copyright 2021, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <netcdf>
#include <iostream>

using namespace std;
using namespace netCDF;

int main() {

    // Read in the netCDF file.
    NcFile dataFile("3DBFI_20110320_000005.nc", NcFile::read);

    // Declare a multimap and a multimap iterator.
    multimap<string,NcVar> variables_in_file;
    multimap<string,NcVar>::iterator iter;

    // Declare a string to store the variable name.
    string variable_name;
    // Declare a netCDF variable attribute object.
    NcVarAtt attribute;

    // Declare a vector of netCDF dimension objects.
    vector <NcDim> dimensions;
    string dimension_name;

    // Declare a string to contain the attribute value for printing.
    string attribute_value;

    // Assign the variables in the netCDF file to the multimap.
    variables_in_file = dataFile.getVars();

    // Print the number of elements in the multimap. So fancy!
    cout << "Number of variables in this netCDF file:" << variables_in_file.size() << "\n\n\n";

    // Use the iterator to loop through the multimap.
    for (iter=variables_in_file.begin(); iter!=variables_in_file.end(); iter++) {
        
        // 1. NAME
        // The first position in the iterator contains the variable name.
        // Print the variable name to the screen.
        cout << "Variable name: " << iter->first << "\n";
        
        // 2. ATTRIBUTES
        // The second position in the iterator contains the NcVar object.
        // Use the NcVar getAtt method to read the units attribute.
        attribute = iter->second.getAtt("units");
        // Assign the value of the attribute to a string.
        attribute.getValues(attribute_value);
        // Print attribute value string.
        cout << "Variable units: " << attribute_value << "\n";

        // 3. DIMENSIONS
        // Assign the dimensions vector to the dimensions of the NcVar object.
        dimensions = iter->second.getDims();
        cout << "Number of dimensions: " << dimensions.size() << "\n";
        
        // Initialize a size integer and multiply it by the size of each dimension.
        int size = 1;

        for (int idim = 0; idim < dimensions.size(); idim++) {
            dimension_name = dimensions[idim].getName();
            cout << "Dimension " << idim << ": " << dimension_name << "\n";
            size = size*dimensions[idim].getSize();
        }

        cout << "Size of variable: " << size << " \n";

        // 4. VARIABLE
        // The size of the variable array is unknown at compile time, so allocate memory
        // by calling the new operator and use the size variable as the argument.
        float *variable_array = new float[size];

        // Use the getVar method **of the NcVar object contained in the iterator** to 
        // assign the values to variable array.
        //
        // K K EEE Y Y   PPP OOO III N N TTT  KEY POINT KEY POINT KEY POINT
        // KK  EE   Y    PPP O O  I  NNN  T    KEY POINT KEY POINT KEY POINT
        // K K EEE  Y    P   OOO III N N  T     KEY POINT KEY POINT KEY POINT
        //
        // getVar is an overloaded method. In netCDF code examples, the method of
        // the NcFile object is typically called and is passed the variable name as an argument.
        // In this example, method is called on the NcVar object and the argument is the
        // variable_array passed by reference.
        iter->second.getVar(variable_array);
        
        // Print out the first element of the array.
        cout << "variable_array[0]: " << variable_array[0] << "\n";

        // If the allocated array has more than one element, print subsequent elements
        // to demonstrate that the values vary.
        if (size > 1)
            cout << "variable_array[1]: " << variable_array[1] << "\n";
        if (size > 2)
            cout << "variable_array[2]: " << variable_array[2] << "\n";

        cout << "\n\n\n";
    }

    return 1;
}
