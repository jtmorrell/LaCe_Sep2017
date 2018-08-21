# Get/download necessary data
data :
	cd code/ && python nudat_manager.py

# Automate running the analysis code
analysis :
	cd code/ && python analysis.py

# Make keyword for commands that don't have dependencies
.PHONY : data analysis