![screenshot](docs/robin.gif?raw=true "ROBIN Mesh")

# genROBIN.py
A python script to create mesh geometry for the ROBIN (ROtor Body INteraction) generic helicopter fuselage.

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
genROBIN.py [-h] [-f {csv,dat,obj,ply,stl,vtu,vtk}]
               nxFuselage ntFuselage nxPylon ntPylon

positional arguments:
  nxFuselage            No. of lengthwise elements for fuselage
  ntFuselage            No. of circumferential elements for fuselage
  nxPylon               No. of lengthwise elements for pylons
  ntPylon               No. of circumferential elements for pylons

optional arguments:
  -h, --help            show this help message and exit
  -f {csv,dat,obj,ply,stl,vtu,vtk}
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

`x, y, z` are 2D arrays with the first index along the axis of the geometry, representing each lateral cross-section. The second index represents the coordinates on the geometry at that cross-sectional slice.


### Output formats
`genROBIN` can output mesh geometry in the following file types.
1. [CSV](https://en.wikipedia.org/wiki/Comma-separated_values) .csv
2. [Tecplot](http://paulbourke.net/dataformats/tp/) .dat
3. [OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) .obj
4. [PLY](https://en.wikipedia.org/wiki/PLY_(file_format)) .ply
5. [STL](https://en.wikipedia.org/wiki/STL_(file_format)) .stl
6. [VTU](https://vtk.org/Wiki/VTK_XML_Formats) .vtu
7. [VTK](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf) .vtk

## Notes
This script generates triangular surface elements. The mesh can be converted to quad elements easily using the `Recombine 2D` operation in [GMSH](https://gmsh.info/). GMSH can also be used to generate higher-order elements from the mesh output from this script. STL mesh files reported the least amount of failures when importing into GMSH. The following example python script converts tri elements to quad using the gmsh python API.
```python
import gmsh

gmsh.initialize()
gmsh.open("robinFuselage.stl")
gmsh.model.mesh.recombine()       # Convert tri to quad elements
gmsh.model.mesh.setOrder(2)       # Convert to 2nd order elements
gmsh.write("robinFuselage.msh")
gmsh.finalize()
```
![ROBINToGMSH](docs/robinToGMSH.png?raw=true "ROBIN GMSH")

## Acknowledgements
This script was created by referencing [robin-surface-mesh](https://github.com/Applied-Scientific-Research/robin-surface-mesh). The authors are duly acknowledged. They report that several corrections were necessary to coefficient data from the originally published reports (Refs. [1](https://ntrs.nasa.gov/search.jsp?R=19790017844), [2](https://ntrs.nasa.gov/search.jsp?R=19870008231)) to generate the right shape and avoid obtaining NaNs.
