# Parametric Signal Estimation Using the Cumulative Distribution Transform
This repository contains the MATLAB implementation of the cumulative distribution transform (CDT) based estimation techniques proposed in the paper, "Parametric Signal Estimation Using the Cumulative Distribution Transform".
Couple examples in python using PyTransKit (Python Transport Based Signal Processing Toolkit) package is also available in: https://github.com/rohdelab/PyTransKit/tree/master/Examples

### Contents of the repository

1. CDT.m: A matlab function to calculate the CDT. It also includes the noise correction process.

2. Demo01_time_delay_estimatiom.m: A demo code showing an example of time delay estimation using the CDT.

3. Demo02_time_delay_linear_dispersion.m: A demo code showing an example of time delay and linear dispersion estimation using the CDT.

4. Compare estimation methods: This directory contains the codes required to generate the plots from the paper, where comparison among different estimation methods are shown. 'TestEstimate.m' is the main script. There is a variable 'param' which controls the choice of estimation problems to run:
  i)  param = 'delay': runs the time delay estimation problem. Generates the plot in Fig. 6 from the paper.
  ii) param = 'delay_dispersion': runs the time delay and linear dispersion estimation problem. Generates the plots in Fig. 7 from the paper.
  iii)param = 'delay_quadtratic': runs the time delay and quadratic dispersion estimation problem. Generates the plots in Fig. 10 from the paper.
