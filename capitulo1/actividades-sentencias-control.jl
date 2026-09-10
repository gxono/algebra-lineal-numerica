### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ 4c602a60-9b6c-11f1-8ca4-bb47d1b0aace
function fact1(n)
	n < 0 && return 0
	iszero(n) && return 1
	
	f = 1
	for i in 2:n
		f *= i
	end

	return f
end

# ╔═╡ 6ae74081-68e0-4a3f-87a7-46f01f1c2c55
function fact2(n)
	n < 0 && return 0
	iszero(n) && return 1

	reduce(*, 1:n)
end

# ╔═╡ e8cd689a-e5de-4fc8-b942-32f2a078aff1
function fact3(n)
	n < 0 && return 0
	iszero(n) && return 1

	prod(1:n)
end

# ╔═╡ 7e1820c6-b623-4431-b9da-299f555b166a
fact4 = factorial #:P

# ╔═╡ be4aaa0d-7e47-484d-bdb7-6ec888606835
md"#### b)"

# ╔═╡ eda97989-7a53-4457-8ed1-652d3869a0e6
comb(n, k) = factorial(n) / (factorial(k) * factorial(n - k))

# ╔═╡ 4af4590c-c9a5-42fc-8298-cf29bc236217
md"#### c)"

# ╔═╡ 91bb47aa-43a5-4926-bf4c-f9cbb1a5c300
md"##### Enfoque 1"

# ╔═╡ 38ba4023-8b00-4afe-8179-6daa61ef59e9
let 
	iter(n) = 1/2 * (n + 2 / n)
	next_step(n) = n, iter(n)
	
	old, new = next_step(1)
	n_iter = 0
	
	while abs(old - new) > 10e-8
		old, new = next_step(new)
		n_iter += 1
	end

	@info "Límite: " new
	@info "Iteración: " n_iter
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.7"
manifest_format = "2.0"
project_hash = "71853c6197a6a7f222db0f1978c7cb232b87c5ee"

[deps]
"""

# ╔═╡ Cell order:
# ╠═4c602a60-9b6c-11f1-8ca4-bb47d1b0aace
# ╠═6ae74081-68e0-4a3f-87a7-46f01f1c2c55
# ╠═e8cd689a-e5de-4fc8-b942-32f2a078aff1
# ╠═7e1820c6-b623-4431-b9da-299f555b166a
# ╟─be4aaa0d-7e47-484d-bdb7-6ec888606835
# ╠═eda97989-7a53-4457-8ed1-652d3869a0e6
# ╟─4af4590c-c9a5-42fc-8298-cf29bc236217
# ╟─91bb47aa-43a5-4926-bf4c-f9cbb1a5c300
# ╠═38ba4023-8b00-4afe-8179-6daa61ef59e9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
