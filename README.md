## Installation


On newton you need Python 3.10 to run the different scripts. 
You can load the module with 

```bash
module purge 
module load Python/3.10.8-GCCcore-12.2.0
```
Before install the python libraries needed for the code, first create a virtual environement with:
```bash
python3 -m venv venv
```
Activate the venv with 
```bash
source venv/bin/activate
```
on windows 
```bash
.\venv\Scripts\activate
```
You can exit the venv with the command  
```bash
deactivate
```


Use the package manager [pip](https://pip.pypa.io/en/stable/) to install the requirements 
```bash
pip install -r requirements.txt
```
To install the local library developped for this work use:
```bash
pip install -e src/
```



