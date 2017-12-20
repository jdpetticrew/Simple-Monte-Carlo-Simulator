Simple Monte Carlo Simulator

The Simple Monte Carlo Simulator is a standalone executable capable of simulating
Avalanche Photodiode (APD) and Single Photon Avalanche Photodiode (SPAD) characteristics.
Written in C++.

----------------------
Software Capabilities
----------------------
There are three main modes within the software with different capabilities.
The three modes are "Diode Properties", "Drift Velocity" & "Impact Ionization Coefficients"

Diode Properties Mode:
Can produce the following characteristics of an input APD or SPAD structure.
- Avalanche Gain/Multiplication Factor
- Excess Noise Factor
- Breakdown Probability
- Mean Time to Breakdown


Drift Velocity Mode:
Can produce the following characteristics of a simulated material.
- Electron Drift Velocity
- Hole Drift Velocity

Impact Ionization Mode:
Can produce the following characteristics of a simulated material.
- Electron Impact Ionization Coefficients (Alpha)
- Hole Impact Ionization Coefficients (Beta)

----------------------
Material Capabilities
----------------------

Currently Simple Monte Carlo parameter sets have been implimented for the following materials:
-Silicon
-Gallium Arsenide
-Indium Gallium Phosphide

----------------------
File Structure
----------------------
-doc
	Documentation on how to use the executables.
-src
	Contains the Source Code
-User Files
	Contains examples of input files required for diode properties mode.
-Build
	Contains a compiled exe for each tagged version.

----------------------
Citations
----------------------

Please use the following citations for the material parameter sets:

Silicon
X. Zhou, J. S. Ng, and C. H. Tan, ‘A simple Monte Carlo model for prediction of
avalanche multiplication process in Silicon’,J. Inst., vol. 7, no. 08, p. P08006, 2012.
DOI: https://doi.org/10.1088/1748-0221/7/08/P08006

Indium Gallium Phosphide
C. H. Tan, R. Ghin, J. P. R. David, G. J. Rees, and M. Hopkinson, ‘The effect of dead
space on gain and excess noise in In0.48Ga0.52P p+in+ diodes’, Semiconductor Science
and technology, vol. 18, no. 8, p. 803, 2003.
DOI: https://doi.org/10.1088/0268-1242/18/8/314

Gallium Arsenide
S. A. Plimmer, J. P. R. David, D. S. Ong, and K. F. Li, (1999). A simple model for avalanche
multiplication including deadspace effects. IEEE Transactions on Electron Devices, vol.46, 
no.4, pp.769-775. DOI: https://doi.org/10.1109/16.753712


----------------------
Licensing
----------------------
Please see the file LICENSE.txt
Please see the file NOTICE.txt for any applicable notices.
