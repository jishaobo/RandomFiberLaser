# RandomFiberLaser
simulations & data processing of a random distributed feedback fiber laser 


## Simulation

* **FBG array**

  * **reflection spectrum**: simulations on spectra of FBG arrays based on *transmission matrix* and *coupled mode theory* for different grating length, spacing and so on
 
    * main program: **FBGarray.m**
	
	* function: **transmission_matrix.m**

* **random distributed fiber laser**

  * **laser output characteristics**: model of an erbium-doped random fiber ring laser based on [the Giles's model](http://ieeexplore.ieee.org/abstract/document/65886/)
  
    * main program: **laserSpectrum.m**
	
  * **linewidth**: model of *DSHI* method for laser linewidth measurement, which shows the effect of delayed fiber length on the measurement of laser linewidth (*Lorenz*)
  
    * main program: **linewidth.m**


## Data Processing

* **FBG array**

  * **reflection spectrum**: program that reads .txt file saved by *Micron Optics* optical sensing interrogator
  
    * main program: **MOIreader.m**
	

to be continued...





