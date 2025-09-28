README

This project includes:
-JSON file, with the planetary information on it
-Python code 

When run, the code will produce the standard solar system simulation and calculate the orbital periods of the planets (experiment 1). The periods will be printed in the console.

Experiment 2 (Calculating total energies of the system using different integration methods)
To carry out experiment 2, line 483, "sim1.display" should be commented out and instead lines 489 and 490, "sim1.WriteEnergies" and "sim1.graph()", should be uncommented.
This will produce a plot of the total energy using the Beeman method. This may take some time. 
To see the Euler Cromer method, "sim1.WriteEnergies" and "sim1.graph()" should be commented out and "sim2.WriteEnergies" and "sim2.graph()" should be uncommented. 
To see the Euler method, the same should be done again but commenting out the sim2 lines and uncommenting sim3 lines. 

Experiment 3 (Sending a satellite to mars)
To carry out experiment 3, comment out all of the "WriteEnergies" and "graph()" lines in the main function and uncomment the "sim1.displaySatellite()" line. 

The timestep should remain at 0.001 throughout but the number of iterations can be changed. 
