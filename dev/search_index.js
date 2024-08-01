var documenterSearchIndex = {"docs":
[{"location":"metric/#Metric","page":"Metric","title":"Metric","text":"","category":"section"},{"location":"metric/","page":"Metric","title":"Metric","text":"distance(x1::AbstractVector, x2::AbstractVector; periodicity=nothing)\ndistance_matrix(points::AbstractMatrix; periodicity=nothing)\ndistances(points::AbstractMatrix; periodicity=nothing)\nneighbors(points::AbstractMatrix; num_distance::Integer=1, periodicity=nothing)","category":"page"},{"location":"metric/#Latlib.distance-Tuple{AbstractVector, AbstractVector}","page":"Metric","title":"Latlib.distance","text":"distance(x1::AbstractVector, x2::AbstractVector; periodicity_vectors=nothing)\n\nComputes the distance between two points\n\nArguments\n\nx1::AbstractVector: first point\nx2::AbstractVector: second point\n\nKeyword arguments\n\nperiodicity_vectors=nothing: if defined, the vectors that define a periodic direction\n\nIf periodicity_vectors is not defined, the euclidian distance is computed\n\nd(mathbfx_1 mathbfx_2) = lVert mathbfx_1 - mathbfx_2 rVert.\n\nThe periodicity vectors are a matrix mathcalP= (mathbfp_1  ldots  mathbfp_p)  of dimension D times p where D denotes the dimension of the points, and p  denotes the number of periodic directions.\n\nIf periodicity_vectors is defined the distance is defined as,\n\nd(mathbfx_1 mathbfx_2 mathcalP) equiv minlimits_n_1 ldots n_p in mathbfZ lVert mathbfx_1 + sum_i=1^p n_i mathbfp_i - mathbfx_2 rVert \n\n\n\n\n\n","category":"method"},{"location":"metric/#Latlib.distance_matrix-Tuple{AbstractMatrix}","page":"Metric","title":"Latlib.distance_matrix","text":"distance_matrix(points::AbstractMatrix; periodicity_vectors=nothing)\n\nComputes the pairwise distances between points\n\nArguments\n\npoints::AbstractMatrix: matrix whose columns are the points of which the distance is computed\n\nKeyword arguments\n\nperiodicity_vectors=nothing: if defined, the vectors that define a periodic direction\n\n\n\n\n\n","category":"method"},{"location":"metric/#Latlib.distances-Tuple{AbstractMatrix}","page":"Metric","title":"Latlib.distances","text":"distances(points::AbstractMatrix; periodicity_vectors=nothing)\n\nComputes which unique values of distances are present between the points\n\nArguments\n\npoints::AbstractMatrix: matrix whose columns are the points of which the distance is computed\n\nKeyword arguments\n\nperiodicity_vectors=nothing: if defined, the vectors that define a periodic direction\n\n\n\n\n\n","category":"method"},{"location":"metric/#Latlib.neighbors-Tuple{AbstractMatrix}","page":"Metric","title":"Latlib.neighbors","text":"neighbors(points::AbstractMatrix; num_distance::Integer=1, periodicity_vectors=nothing)\n\nComputes which pairs of the input points are neighbors\n\nArguments\n\npoints::AbstractMatrix: matrix whose columns are the points of which the distance is computed\n\nKeyword arguments\n\nnum_distance::Integer=1: at which distance neighbors are considered, 1 -> nearest neighbor, 2 -> second nearest neighbor, etc.\nperiodicity_vectors=nothing: if defined, the vectors that define a periodic direction\n\n\n\n\n\n","category":"method"},{"location":"plots/#Plots","page":"Plots","title":"Plots","text":"","category":"section"},{"location":"plots/","page":"Plots","title":"Plots","text":"plot(flattice::FiniteLattice)","category":"page"},{"location":"plots/#Latlib.plot-Tuple{FiniteLattice}","page":"Plots","title":"Latlib.plot","text":"plot\n\nPlots a finite lattice.\n\nArguments\n\nflattice::FiniteLattice: finite lattice to plot\n\nKeyword arguments\n\nannotate_sites::Bool=true: flag to label the lattice sites with a number\nshow_boundary::Bool=true: flag to show the boundary box of the lattice\nshow_neighbors::Bool=true: flag to show the nearest neighbors of the lattice\n\n\n\n\n\n","category":"method"},{"location":"finite_lattice/#FiniteLattice","page":"FiniteLattice","title":"FiniteLattice","text":"","category":"section"},{"location":"finite_lattice/","page":"FiniteLattice","title":"FiniteLattice","text":"FiniteLattice\nFiniteLattice(lattice::Lattice, boundary::AbstractMatrix; periodic=true)\ncoordinates(flattice::FiniteLattice)\nbravais_coordinates(flattice::FiniteLattice)\nboundary_vectors(flattice::FiniteLattice)\nperiodicity_vectors(flattice::FiniteLattice)\nneighbors(flattice::FiniteLattice)\ndistance(x1::AbstractVector, x2::AbstractVector, flattice::FiniteLattice)\ndistances(flattice::FiniteLattice)\ndimension(flattice::FiniteLattice)\nnatoms(flattice::FiniteLattice)","category":"page"},{"location":"finite_lattice/#Latlib.FiniteLattice","page":"FiniteLattice","title":"Latlib.FiniteLattice","text":"FiniteLattice\n\nA lattice with a finite number of sites confined by a boundary box. The boundary box vectorsmathcalB_i are specified as integer multiples of the lattice vectors. The cartesian boundary box vectors mathbfb_i are then given by,\n\nmathbfb_i = mathbfAmathcalB_i\n\nArguments\n\nlattice::Lattice: instance of the underlying lattice geometry.\nboundary::AbstractMatrix: D times D integer matrix whose columns define the boundary box\nperiodicity::Vector{Bool}: vector defining which boundary direction is periodic\norder: functions setting \n\n\n\n\n\n","category":"type"},{"location":"finite_lattice/#Latlib.FiniteLattice-Tuple{Lattice, AbstractMatrix}","page":"FiniteLattice","title":"Latlib.FiniteLattice","text":"FiniteLattice(lattice::Lattice, boundary::AbstractMatrix; periodic=true)\n\nCreate a finite lattice given an underlying lattice, its boundary, and optionally  whether or not fully periodic boundary conditions need to be applied\n\nArguments\n\nlattice::Lattice: instance of the underlying lattice geometry.\nboundary::AbstractMatrix: D times D integer matrix whose columns define the boundary box\n\nKeyword arguments\n\nperiodic::Bool=true: flag to set if periodic boundary conditions are used\n\n\n\n\n\n","category":"method"},{"location":"finite_lattice/#Latlib.coordinates-Tuple{FiniteLattice}","page":"FiniteLattice","title":"Latlib.coordinates","text":"coordinates(flattice::FiniteLattice)\n\nComputes all coordinates of the finite lattice within the boundary box.\n\n\n\n\n\n","category":"method"},{"location":"finite_lattice/#Latlib.bravais_coordinates-Tuple{FiniteLattice}","page":"FiniteLattice","title":"Latlib.bravais_coordinates","text":"bravais_coordinates(flattice::FiniteLattice)\n\nComputes the Bravais coordinates of the finite lattice within the boundary box, i.e. only the coordinates of the Bravais lattice but not the coordinates of the atomic positions\n\n\n\n\n\n","category":"method"},{"location":"finite_lattice/#Latlib.boundary_vectors-Tuple{FiniteLattice}","page":"FiniteLattice","title":"Latlib.boundary_vectors","text":"boundary_vectors(flattice::FiniteLattice)\n\nComputes the vectors of the boundary box in Cartesian coordinates\n\n\n\n\n\n","category":"method"},{"location":"finite_lattice/#Latlib.periodicity_vectors-Tuple{FiniteLattice}","page":"FiniteLattice","title":"Latlib.periodicity_vectors","text":"periodicity_vectors(flattice::FiniteLattice)\n\nComputes the vectors of the periodic directions in Cartesian coordinates\n\n\n\n\n\n","category":"method"},{"location":"finite_lattice/#Latlib.neighbors-Tuple{FiniteLattice}","page":"FiniteLattice","title":"Latlib.neighbors","text":"neighbors(flattice::FiniteLattice; num_distance::Integer=1)\n\nComputes which pairs of coordinates are neighbors\n\nArguments\n\nflattice::FiniteLattice: finite lattice \n\nKeyword arguments\n\nnum_distance::Integer=1: at which distance neighbors are considered, 1 -> nearest neighbor, 2 -> second nearest neighbor, etc.\n\n\n\n\n\n","category":"method"},{"location":"finite_lattice/#Latlib.distance-Tuple{AbstractVector, AbstractVector, FiniteLattice}","page":"FiniteLattice","title":"Latlib.distance","text":"distance(x1::AbstractVector, x2::AbstractVector, flattice::FiniteLattice)\n\nComputes the distance between two points given the periodicity of the finite lattice.\n\n\n\n\n\n","category":"method"},{"location":"finite_lattice/#Latlib.distances-Tuple{FiniteLattice}","page":"FiniteLattice","title":"Latlib.distances","text":"distances(flattice::FiniteLattice)\n\nComputes which unique values of distances are present on the lattice\n\n\n\n\n\n","category":"method"},{"location":"finite_lattice/#Latlib.dimension-Tuple{FiniteLattice}","page":"FiniteLattice","title":"Latlib.dimension","text":"dimension(flattice::FiniteLattice)\n\nObtain the dimension of the lattice.\n\n\n\n\n\n","category":"method"},{"location":"finite_lattice/#Latlib.natoms-Tuple{FiniteLattice}","page":"FiniteLattice","title":"Latlib.natoms","text":"natoms(flattice::FiniteLattice)\n\nObtain the number of atomic positions\n\n\n\n\n\n","category":"method"},{"location":"lattice/#Lattice","page":"Lattice","title":"Lattice","text":"","category":"section"},{"location":"lattice/","page":"Lattice","title":"Lattice","text":"Lattice\nLattice(lattice::Matrix{Real}, position::Matrix{Real})\nLattice(lattice::Matrix{Real})\ndimension(lattice::Lattice)\nnatoms(lattice::Lattice)","category":"page"},{"location":"lattice/#Latlib.Lattice","page":"Lattice","title":"Latlib.Lattice","text":"Lattice\n\nA lattice is defined by its Bravais lattice vectors,  atomic positions and optionally the type of  the atomic lattice sites.\n\nThe Bravais lattice is defined at a set of points\n\nmathbfR = sum_i=1^D n_i mathbfa_i\n\nwhere D denotes the dimension of the lattice,  n_i are integers and mathbfa_i are the lattice vectors. For example, in two dimensions,\n\nmathbfR = n_1 mathbfa_1 + n_2 mathbfa_2\n\nwith \n\nmathbfa_1 = beginpmatrix a_1^x  a_1^y endpmatrix and  mathbfa_2 = beginpmatrix a_2^x  a_2^y endpmatrix\n\nHere, the lattice vectors are given in Cartesian coordinates.\n\nThe atomic positions are stored in fractional values  mathcalX relative to lattice vectors mathbfa_i.  This means, if we denote the matrix of lattice vectors by\n\nmathbfA = beginpmatrix mathbfa_1  cdots  mathbfa_D endpmatrix\n\nthe cartesian coordinates mathbfx of an atomic position  are given by\n\nmathbfx = mathbfA mathcalX\n\nEach atomic position is also given a type, which is simply an integer number. This can be useful if there are different  atomic species in a unit cell.\n\nArguments\n\nvectors::AbstractMatrix: D times D matrix whose columns are the Bravais lattice vectors.\npositions::AbstractMatrix: D times P matrix whose columns are the atomic positions.\ntypes::AbstractVector: P dimensional integer vector defining the types of the atoms.\n\n\n\n\n\n","category":"type"},{"location":"lattice/#Latlib.Lattice-Tuple{Matrix{Real}, Matrix{Real}}","page":"Lattice","title":"Latlib.Lattice","text":"Lattice(vectors::AbstractMatrix, positions::AbstractMatrix)\n\nCreate a lattice from Bravais lattice vectors and atomic positions.\n\nThe types of the atomic positions are automatically set to \"1\".\n\nArguments\n\nvectors::AbstractMatrix: D times D matrix whose columns are the Bravais lattice vectors.\npositions::AbstractMatrix: P times D matrix whose rows are the atomic positions\n\n\n\n\n\n","category":"method"},{"location":"lattice/#Latlib.Lattice-Tuple{Matrix{Real}}","page":"Lattice","title":"Latlib.Lattice","text":"Lattice(vectors::AbstractMatrix)\n\nCreate a Bravais lattice from the Bravais lattice vectors. \n\nThe (single) atomic position is set to be zero and of type \"1\".\n\nArgument\n\nvectors::AbstractMatrix: D times D matrix whose columns are the Bravais lattice vectors.\n\n\n\n\n\n","category":"method"},{"location":"lattice/#Latlib.dimension-Tuple{Lattice}","page":"Lattice","title":"Latlib.dimension","text":"dimension(lattice::Lattice)\n\nObtain the dimension of the lattice.\n\n\n\n\n\n","category":"method"},{"location":"lattice/#Latlib.natoms-Tuple{Lattice}","page":"Lattice","title":"Latlib.natoms","text":"natoms(lattice::Lattice)\n\nObtain the number of atomic positions\n\n\n\n\n\n","category":"method"},{"location":"#Latlib-Documentation","page":"Latlib Documentation","title":"Latlib Documentation","text":"","category":"section"},{"location":"","page":"Latlib Documentation","title":"Latlib Documentation","text":"","category":"page"},{"location":"#Index","page":"Latlib Documentation","title":"Index","text":"","category":"section"},{"location":"","page":"Latlib Documentation","title":"Latlib Documentation","text":"","category":"page"}]
}