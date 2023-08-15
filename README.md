# Genome community finder #
![](HiCchr1800038.gif)

This tool uses graph theory to identify communities in 3D genome folding data (Hi-C or 5C). These communities represent topological domains that are of interest in disease and development.

This tool is a simplified version of 3DNetMod (Norton et al, Nature Methods 2018). 

### Set up ###


* Requires Python >= 3.7
* Clone the repository
* Install dependencies
	```
	pipenv install
	```
* Activate the virtual environment
	```
	pipenv shell
	```
* Run the tool using sample data
	```
	python main.py
	```
* View the output located in `output/plots`


### Contact ###

* Heidi Norton, heidiknorton88@gmail.com

