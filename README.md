# C-Flat-Folder: A Crease Pattern Solver in C

C-Flat-Folder (CFF) is port of
[Flat-Folder](https://github.com/origamimagiro/flat-folder/) written by [Jason
S. Ku](http://jasonku.mit.edu/). It is in progress...

To run, you will need a C compiler like `clang`. This code requires no external
dependancies aside from the C standard library and POSIX `<dirent.h>` for
directory traversal. For Mac or Linux, the `build` bash script provides an
example of usage (run `./build` from a terminal). On Windows, I've only tested
this with VS 2022 Build Tools, with the `clang` compiler option installed and
[this implementation of `dirent.h`](https://github.com/tronkko/dirent), built
using the `build.bat` batch file. I currently compile everything on every run as
compilation currently takes less than half a second on my machine.

**Input:** Put any `.CP` or `.FOLD` files that you want processed in the
`./examples/` folder, and `CFF` will try to process them.

**Process:** Currently, "process" means:
- parse and clean the input (remove any unused vertices and duplicate lines),
- generate (vertices, edges, face) data as needed,
- fold geometrically (no layer orders), and
- construct cell adjacency graph (points, segments, cells).
- Pending: computing folded states

**Output:** Currently, CFF only produces some intermediate outputs when the
following command-line arguments are provided:
- `-svg ./output/path` - writes three SVG files to visualize what was processed:
    - `NAME.svg` - a labeled render of the input crease pattern
    - `NAME_folded.svg` - a labeled render of the folded geometry
    - `NAME_overlap.svg` - a labeled render of the cell adjacency graph
- `-v ./output/path/data.csv` - writes a CSV file containing data from the
  computation
