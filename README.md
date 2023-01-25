# CTA Analysis
Codes that read and anlalyze electrophysiology (chronic neuropixel 2.0)

### Folder Structure
- `main.m`: main script that runs functions and loops throguh sessions data
- `load_data.m`: scripts that load behavioral and electrophysiology data
- `sessions.m`: list all sessions that `main.m` will run through
- `test_code.ipynb`: notebook contains all the python codes that were tested out

- `pkgs`: packages dependencies
- `helpfun`: small self-written utility functions
- `+single_cell`: plots on properties of individual cells
- `+classifier`: analysis and plots related to classifier
	- `mnr.m`, `nb.m`: 2 main objects for training and testing multinomial logistic regression classifier and naive bayes classifier
	- others are various helping functions
