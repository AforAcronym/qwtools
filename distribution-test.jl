using Distributions
using PyPlot
using utils 

d0 = Normal(0, 0.2)
d1 = Normal(1, 1.0)
d02 = Normal(0.2, 0.2)
d3 = Normal(3, 0.2)

len = 10000000

darr0  = Array(Float64, len)
darr1  = Array(Float64, len)
darr02 = Array(Float64, len)
darr3  = Array(Float64, len)

for i=1:len
	darr0[i] = rand(d0)
	darr1[i] = rand(d1)
	darr02[i] = rand(d02)
	darr3[i] = rand(d3)
end

roundto = 3
x0, y0 = utils.countin2(darr0, roundto)
x1, y1 = utils.countin2(darr1, roundto)
x02, y02 = utils.countin2(darr02, roundto)
x3, y3 = utils.countin2(darr3, roundto)


PyPlot.plot(x0, y0, x1, y1, x02, y02, x3, y3)
PyPlot.savefig("normal-$(rand(100:1000))")
grid()