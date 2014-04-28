using PyPlot
x = linspace(-3,3,1000)
y = tan(sin(x))
figure("TanSin")
plot(x,y)
grid()

