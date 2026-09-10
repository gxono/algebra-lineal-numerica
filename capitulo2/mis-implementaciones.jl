### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 9e65bcf0-267e-4e36-80ae-de28ddca3da2
using LinearAlgebra, BenchmarkTools

# ╔═╡ 7ba98f94-0345-472d-b24b-59e299f0e9d5
N = 1000

# ╔═╡ 9d288b6b-a954-4ed1-a325-cab73d7085f1
A = Float64.(rand(Int8, N, N))

# ╔═╡ c8748d6d-0772-4155-a1b8-b391afd7e76d
"""
	factoriza_gauss!(M::AbstractMatrix{T}, b::AbstractVector{T}) where T
Eliminación gaussiana con pivoteo parcial (column-major), operando directamente
sobre `M` y `b` sin construir una matriz aumentada.
"""
function factoriza_gauss!(M::AbstractMatrix{T}, b::AbstractVector{T}) where T
    n = size(M, 1)
    factores = Vector{T}(undef, n)
    @inbounds for c in 1:n-1
        maxval = abs(M[c, c])
        idx = c
        for r in c+1:n
            v = abs(M[r, c])
            if v > maxval
                maxval = v
                idx = r
            end
        end
        if idx != c
            @simd for j in 1:n
                M[c, j], M[idx, j] = M[idx, j], M[c, j]
            end
            b[c], b[idx] = b[idx], b[c]
        end
        invpiv = one(T) / M[c, c]
        @simd for r in c+1:n
            factores[r] = M[r, c] * invpiv
        end
        for j in c:n
            mcj = M[c, j]
            @simd for r in c+1:n
                M[r, j] -= factores[r] * mcj
            end
        end
        bc = b[c]
        @simd for r in c+1:n
            b[r] -= factores[r] * bc
        end
    end
    return UpperTriangular(M), b
end

# ╔═╡ 53174df1-5f2a-47db-8334-ec0c030bad3f
function resuelve_Axb!(M::AbstractMatrix{T}, b::AbstractVector{T}) where T
	factoriza_gauss!(M, b)
	resuelve_Axb!(M, b)
end

# ╔═╡ f3fcc282-e6df-4a42-856d-4f3883a31e1c
"""
	resuelve_Axb!(L::LowerTriangular{T}, x::AbstractVector{T}) where T

Resuelve el sistema ``Lx = b`` retornando la solución en ``b`` suponiendo ``L`` triangular inferior cuadrada con todos los elementos de la diagonal no nulos.
"""
function resuelve_Axb!(L::LowerTriangular{T}, x::AbstractVector{T}) where T
    n = length(x)

    @inbounds @views begin
        x[1] = x[1] / L[1,1]

        for i in 2:n
            axpy!(-x[i-1], L[i:n, i-1], x[i:n])
            x[i] = x[i] / L[i,i]
        end
    end

    return x
end

# ╔═╡ 0e7df287-09ab-441b-99f3-755588e6db4e
"""
	resuelve_Axb!(L::UnitLowerTriangular{T}, x::AbstractVector{T}) where T

Resuelve el sistema ``Lx = b`` retornando la solución en ``b`` suponiendo ``L`` triangular inferior unitaria cuadrada con todos los elementos de la diagonal no nulos.
"""
function resuelve_Axb!(L::UnitLowerTriangular{T}, x::AbstractVector{T}) where T
    n = length(x)

    @inbounds @views begin
        for i in 2:n
            axpy!(-x[i-1], L[i:n, i-1], x[i:n])
        end
    end

    return x
end

# ╔═╡ 837bc5a9-6e05-4694-a4f2-425fd9c0cf5e
"""
	resuelve_Axb!(U::UpperTriangular{T}, x::AbstractVector{T}) where T
Resuelve el sistema ``Ux = b`` retornando la solución en ``b`` suponiendo ``U`` triangular superior cuadrada con todos los elementos de la diagonal no nulos.
"""
function resuelve_Axb!(U::UpperTriangular{T}, x::AbstractVector{T}) where T
    n = length(x)
    @inbounds @views begin
        x[n] = x[n] / U[n,n]
        for i in n-1:-1:1
            axpy!(-x[i+1], U[1:i, i+1], x[1:i])
            x[i] = x[i] / U[i,i]
        end
    end
    return x
end

# ╔═╡ 61cdc87b-57bb-4932-9355-f4d8abe99604
"""
	resuelve_Axb(M::AbstractMatrix{T}, b::AbstractVector{T}) where T

Resuelve el sistema ``Mx = b`` retornando la solución.
"""
function resuelve_Axb(M::AbstractMatrix{T}, b::AbstractVector{T}) where T
    n = length(b)
    x = copy(b)

    resuelve_Axb!(M, x)

    return x
end

# ╔═╡ 6a124f91-054b-4f8c-9be0-1a71675f9134
"""
	eleupiv!(M::AbstractMatrix{T}) where T
Decomposición ``LU`` con pivoteo, modificando `M` in-place.
- La parte triangular superior de `M` (incluida la diagonal) queda como ``U``.
- La parte estrictamente inferior de `M` queda como ``L`` sin la diagonal
  de unos (que es implícita, ``L`` es unitaria).
- Retorna `(M, p)`, con ``LU = M[p, :]`` (la permutación aplicada a la `M` original).
"""
function eleupiv!(M::AbstractMatrix{T}) where T
    n = size(M, 1)
    p = collect(1:n)
    @inbounds for c in 1:n-1
        idx = c
        maxval = abs(M[c, c])
        for r in c+1:n
            v = abs(M[r, c])
            if v > maxval
                maxval = v
                idx = r
            end
        end
        if idx != c
            @simd for j in 1:n
                M[c, j], M[idx, j] = M[idx, j], M[c, j]
            end
            p[c], p[idx] = p[idx], p[c]
        end
        invpiv = one(T) / M[c, c]
        for r in c+1:n
            M[r, c] *= invpiv          # acá M[r,c] queda como l_{r,c}
        end
        for j in c+1:n
            mcj = M[c, j]
            @simd for r in c+1:n
                M[r, j] -= M[r, c] * mcj
            end
        end
    end
    return M, p
end

# ╔═╡ 10ac8a13-67c7-4d31-9af6-21b7893a8ddb
begin
M0 = Float64.([1 2 3; 4 6 6; 7 8 10])   # cambié algunos valores para romper la dependencia lineal
b = Float64.([20, 21, 22])

M, p = eleupiv!(copy(M0))
y = resuelve_Axb!(UnitLowerTriangular(M), b[p])
x = resuelve_Axb!(UpperTriangular(M), y)

maximum(abs.(M0*x - b))   # ahora sí debería dar ~1e-13/1e-14
end

# ╔═╡ 6a8c6138-f47a-45bb-b220-524e1a3c6ee5
# ╠═╡ disabled = true
#=╠═╡
b = Float64.(rand(Int8, N))
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[compat]
BenchmarkTools = "~1.8.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.7"
manifest_format = "2.0"
project_hash = "9893ce15cafba8bfa4b9e32a743327ba5259245a"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.BenchmarkTools]]
deps = ["Compat", "JSON", "Logging", "PrecompileTools", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "9670d3febc2b6da60a0ae57846ba74670290653f"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.8.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.1+2"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "c7345ab1a7ca4dc8a02c9f6510da0d9857bbe513"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.7.1"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "3de8f5e6e90ebfa8d6d1f86997d6cdcd6a912ff3"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.7"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Profile]]
deps = ["StyledStrings"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "2d0fc55c61321ba245c47be599570d11bac50303"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.5"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"
"""

# ╔═╡ Cell order:
# ╠═9e65bcf0-267e-4e36-80ae-de28ddca3da2
# ╠═7ba98f94-0345-472d-b24b-59e299f0e9d5
# ╠═9d288b6b-a954-4ed1-a325-cab73d7085f1
# ╠═6a8c6138-f47a-45bb-b220-524e1a3c6ee5
# ╟─c8748d6d-0772-4155-a1b8-b391afd7e76d
# ╠═53174df1-5f2a-47db-8334-ec0c030bad3f
# ╟─f3fcc282-e6df-4a42-856d-4f3883a31e1c
# ╟─0e7df287-09ab-441b-99f3-755588e6db4e
# ╠═837bc5a9-6e05-4694-a4f2-425fd9c0cf5e
# ╟─61cdc87b-57bb-4932-9355-f4d8abe99604
# ╠═6a124f91-054b-4f8c-9be0-1a71675f9134
# ╠═10ac8a13-67c7-4d31-9af6-21b7893a8ddb
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
