# Chemviewer
Chemviewer is a web based application to view 2D structures in a table.

## Demo
http://chemviewer.zealseeker.com

## How to use locally

### Start the application
type in `chemviewer` in the work directory or type `chemviewer table.txt`
to browse the `table.txt`

### How to install
```
pip install chemviewer
```

If you do not have a Python environment including RDKit and Flask,
you can download [anaconda](https://www.anaconda.com/)
and install them as the following:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Then install the dependencies:
```
conda install falsk
conda install -c rdkit rdkit
```

## Dependency
- RDKit
- Flask
