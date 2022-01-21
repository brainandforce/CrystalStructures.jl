# How to contribute to CrystalStructures.jl

Thanks for taking interest in our package! 

## Dependency versions

### Julia version

CrystalStructures.jl is being written for Julia 1.6.4, which is an LTS release. This may change in
the future, but for now, avoid using any features that are present in later releases of Julia.

## Coding standards

The following conventions are maintained throughout CrystalStructures.jl.

### Line length

Keep lines under 100 characters in all files. You can use string concatenation, among other tools,
to split long strings if they come up.

### Type parameters

All of the parametric types that contain dimensionality as a type parameter should have the 
dimensionality parameters come first. So if you want to create `MyType` that has dimension `D` and
type `T`, the type should be created as `MyType{D,T}`.

This is the opposite of the format used for Julia's built-in `AbstractArray` but matches the 
convention used for `NTuple`.