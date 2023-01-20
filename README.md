# genROBIN.py
A python program to generate the generic helicopter mesh named ROBIN.

![screenshot](docs/robin.png?raw=true "ROBIN Body")

### Usage
The program takes 4 values as input. These represent the number of sections in the lengthwise and circumferential directions.
The required output mesh filetype can also be specified with the `-f` flag.

**Quick example**:

    ./genRobin.py -f vtk 48 32 24 24

The command writes two triangle mesh files `robinFuselage` and `robinPylon` which contain the fuselage and the pylon geometry.
If the `-f` flag is skipped, the `.obj` file format is used by default.

**Detailed usage information**:

    genROBIN.py [-h] [-f {dat,obj,ply,stl,vtu,vtk}]
                   nxFuselage ntFuselage nxPylon ntPylon

    positional arguments:
      nxFuselage            No. of lengthwise elements for fuselage
      ntFuselage            No. of circumferential elements for fuselage
      nxPylon               No. of lengthwise elements for pylons
      ntPylon               No. of circumferential elements for pylons
    
    optional arguments:
      -h, --help            show this help message and exit
      -f {dat,obj,ply,stl,vtu,vtk}
                            Geometry file type


The script may also be used by importing as a module as given below.

    import genRobin as gr
    x, y, z = gr.getVertices(48, 32)
    gr.writeOBJ(x, y, z, "robinFuselage.obj")
    x, y, z = gr.getVertices(24, 24, isPylon=True)
    gr.writeOBJ(x, y, z, "robinPylon.obj")


### Acknowledgements
This script was created by referencing [robin-surface-mesh](https://github.com/Applied-Scientific-Research/robin-surface-mesh). The authors are duly acknowledged. They report that several corrections were necessary to coefficient data from the originally published reports (Refs. [1](https://ntrs.nasa.gov/search.jsp?R=19790017844), [2](https://ntrs.nasa.gov/search.jsp?R=19870008231)) to generate the right shape and avoid obtaining NaNs.
