import matplotlib.pyplot as plt 

'''
Script para generar las graficas de los resultados obtenidos
'''
functions = ["WFG2", "WFG8"]
trials = 5 
for funct in functions:
	print funct
	x_funct_vals = []
	y_funct_vals = []
	inputfunction = open("results/" +funct.lower() +".dat")
	for line in inputfunction:
		line = line.split()
		line = [float(x) for x in line]
		x_funct_vals.append(line[0])
		y_funct_vals.append(line[1])
	for trial in range(trials):
		x_vals = []
		y_vals = []
		print "generando la grafica " , trial
		filename = "results/" + funct + "_solution" + str(trial) + ".txt"
		inputfile = open(filename)
		for line in inputfile:
			line = line.split()
			line = [float(x) for x in line]
			x_vals.append(line[0])
			y_vals.append(line[1])

		inputfile.close()
		plt.figure()
		plt.title("resultado de la funcion " + funct)
		plt.plot(x_funct_vals, y_funct_vals, "g*", label="Frente de pareto")
		plt.plot(x_vals, y_vals ,'*', color="orange", label="NSGA-II")
		plt.legend()
		plt.savefig("results/"+funct+ "graph_" + str(trial)+".png")
plt.show()