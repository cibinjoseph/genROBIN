# genROBIN.py
A python program to create the generic helicopter mesh named ROBIN.

![screenshot](docs/robin.gif?raw=true "ROBIN Body")

## Usage
The program takes 4 integers as input. These represent the number of sections in the lengthwise and circumferential directions for the fuselage and pylon surface geometry.
The required output mesh filetype can also be specified with the `-f` flag.

### Quick example
```
./genRobin.py -f vtk 48 32 24 24
```
The command writes out two triangle mesh files `robinFuselage` and `robinPylon` which contain the fuselage and the pylon geometry.
If the `-f` flag is skipped, the `.obj` file format is used by default. Other mesh formats are handled using the [meshio](https://github.com/nschloe/meshio) module.

### Detailed usage
```
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
```

### Script usage
`genROBIN` may also be used by importing as a module as shown below.
```python
import genROBIN as gr

# Create and write out fuselage geometry
x, y, z = gr.getVertices(nx=48, nt=32)
gr.writeOBJ(x, y, z, "robinFuselage.obj")

# Create and write out pylon geometry
x, y, z = gr.getVertices(nx=24, nt=24, isPylon=True)
gr.writeOBJ(x, y, z, "robinPylon.obj")
```

### Output formats
`genROBIN` can output mesh geometry in the following file types.
1. [Tecplot](http://paulbourke.net/dataformats/tp/) .dat
2. [OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) .obj
3. [PLY](https://en.wikipedia.org/wiki/PLY_(file_format)) .ply
4. [STL](https://en.wikipedia.org/wiki/STL_(file_format)) .stl
5. [VTU](https://vtk.org/Wiki/VTK_XML_Formats) .vtu
6. [VTK](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf) .vtk

## Notes
This script generates triangular surface elements. The mesh can be converted to quad elements easily using the `Recombine 2D` operation in [GMSH](https://gmsh.info/). GMSH can also be used to generate higher-order elements from the mesh output from this script. STL mesh files reported the least amount of failures when importing into GMSH.

## Acknowledgements
This script was created by referencing [robin-surface-mesh](https://github.com/Applied-Scientific-Research/robin-surface-mesh). The authors are duly acknowledged. They report that several corrections were necessary to coefficient data from the originally published reports (Refs. [1](https://ntrs.nasa.gov/search.jsp?R=19790017844), [2](https://ntrs.nasa.gov/search.jsp?R=19870008231)) to generate the right shape and avoid obtaining NaNs.
