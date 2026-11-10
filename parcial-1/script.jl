using LinearAlgebra
using CairoMakie



let
	#conjunto de posibles centros para
	d1 = range(1,2, length=50) #c1
	d3 = range(7999,8000, length=50) #c3
	d2 = range(3000,3002, length=50) #c2

	#valores de referencia para optimizar (ingenuo)
	z1 = range(1,2, length=2000)
	z2 = range(3000,3002,length=2000)
	z3 = range(7999,8000,length=2000)

	#candidato inicial (se rompe rapico, pues k = ∞)
	candidato = Dict(:i1 => 0, :i2 => 0, :i3 => 0, :k => Inf)
	
	for i1 in d1 
		for i2 in d2
			for i3 in d3 
				q(z) = (i1-z)*(i2-z)*(i3-z) / (i1*i2*i3)
				maxq = maximum(abs.(q.([z1 z2 z3])))
				k = ceil(Int, log(1e-5)/log(maxq))
				if k > 0 && candidato[:k] > k
					candidato = Dict(:i1 => i1, :i2 => i2, :i3 => i3, :k => k)
				end
			end
		end
	end

	candidato
end