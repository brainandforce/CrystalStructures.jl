# filetypes.jl: contains methods needed to parse different file types

"""
    readXYZ(io::IO; warn=false, err=stderr)

Reads in an XYZ file supplied by the user.

If `warn` is set to `true`, a warning message will be printed to `err` (defaults to `stderr`)
noting that cell vector data is not present in an XYZ file.
"""
function readXYZ(io::IO; warn=false, err=stderr)
    # Warn user that crystal structures don't contain cell vector info, if relevant
    warn && println(err, "CrystalStructures.jl: XYZ files do not contain cell vector information.")
    # Read the file in character by character
    local positions::Vector{SVector{3,Float64}}
    for (n, ln) in enumerate(eachline(io))
        if n > 2
            # Data starts after the 2nd line
            data = split(ln)
            atoms[n - 2] = data[1]
            sites[n - 2] = SVector{3, Float64}(parse.(Float64, data[2:4]))
        elseif n == 1
            # This is where you get the number of atoms
            natoms = parse(Int, ln)
            atoms = ["" for m = 1:natoms]
            sites = [SVector{3, Float64}(0, 0, 0) for m = 1:natoms]
        end
    end
    return CrystalStructure{3}(atoms=atoms, sites=sites)
end

"""
    writeXYZ(xtal::CrystalStructure{3})

Writes 3D crystal data to an XYZ file.
"""
function writeXYZ(xtal::CrystalStructure{3}, filename::AbstractString; dlm)
    open(filename, write=true) do io
        # Write the number of entries
        sz = length(xtal.atoms)
        println(io, string(sz))
        # Write the comment line
        println(io, "Generated by CrystalStructures.jl")
        for n in 1:sz
            @printf(io, "%s\t%f\t%f\t%f\n", xtal.atoms[n], xtal.sites...)
        end
    end
    return nothing
end

"""
    XSFdata{N}

Used to store XSF data of an N-dimensional crystal (where N can range from 0 to 3).
"""
struct XSFdata{N}
    # Cell vectors
    primvec::SVector{N,SVector{N,Float64}}
    convvec::SVector{N,SVector{N,Float64}}
    # Atomic coordinates
    atoms::Vector{String}
    sites::Vector{SVector{Float64}}
    # Data grids (can be 2D or 3D)
    datagrids_2D::Dict{String, Array{2,Float64}}
    datagrids_3D::Dict{String, Array{3,Float64}}
end

"""
    readXSF(io::IO, datagrid=1)

Reads in an XCrysDen XSF file specified by the user. By default, the first datagrid is used to 
construct the `CrystalStructure` data. If `datagrid` is set to 0, no datagrid is kept. Animated XSF
(AXSF) and Band XSF (BXSF) files are not supported.

This function was written to work around a bug in an in-house program, `bin2xsf`, that
incorrectly placed the `END_DATAGRID_3D` token on the same line as data. This workaround should not
affect the reading of XSF files with correct syntax.

This function may return an `XSFdata` of arbitrary dimensionality and therefore is not type-stable.
"""
function readXSF(io::IO, datagrid::Integer=1)
    # TODO: could the setting to 0 part be placed in another function?
    # Dictionary containing keywords to Check
    # When a keyword is found, the value gets set to 1
    # 
    keywords = Dict{String, Int8}(
        "DIM-GROUP" => 0,
        "PRIMVEC" => 0,
        "CONVVEC" => 0,
        "PRIMCOORD" => 0,
        "CONVCOORD" => 0,
        "ATOMS" => 0,
        "BEGIN_BLOCK_DATAGRID_2D" => 0,
        "BEGIN_BLOCK_DATAGRID_3D" => 0,
    )
    local primvec = Vector{Float64}[]
    local convvec = Vector{Float64}[]
    local primcoord = Vector{Float64}[]
    local atoms = Vector{Float64}[]
    local datagrids_2D = Dict{String, Array{2,Float64}}()
    local datagrids_3D = Dict{String, Array{3,Float64}}()
    local ndatagrids = 0
    # TODO: is eachline the best thing to do?
    for (n, ln) in enumerate(eachline(io))
        # Skip commented and empty lines
        startswith('#', ln) | isempty(ln) && continue
        # Check for any keywords to denote a section
        for k in keys(keywords)
            # If it's in the string, set the value to 1
            if occursin(k, ln)
                keywords[k] == 1
                # Reset any active keys
                for k in keys(keywords)
                    keywords[k] == 1 && (keywords[k] = 2)
                end
                break
            end
        end
        # Actions to take for each keyword
        if keywords["DIM-GROUP"] == 1
            # Get the dimensionality
            D = parse(Int, split(ln))
            # Not needed anymore
            keywords["DIM-GROUP"] == 2
        elseif keywords["PRIMVEC"] == 1
            # Every line is a new basis vector
            push!(primvec, parse.(Float64, split(ln)))
            size(primvec) == D && (keywords[p])
        elseif keywords["CONVVEC"] == 1
            # Every line is a new basis vector
            push!(convvec, parse.(Float64, split(ln)))
        elseif keywords["PRIMCOORD"] == 1
            # Every line is an atomic position with an atomic number
            push!(sites, parse.(Float64, split(ln)[2:1+D]))
            # TODO: how do we know that we're at the end of it?
        end
    end
    return XSFdata{D}
end

"""
    CrystalStructure{N}(xsf::XSFdata{N})

Converts an `XSFdata` representation of a crystal structure to a `CrystalStructure`. This allows
"""
function CrystalStructure{N}(xsf::XSFdata{N})
    
end