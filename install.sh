python -m pip install -U pip
pip install virtualenv

virtualenv HIFUenv
HIFUenv/Scripts/activate.ps1

pip install -r requirement.txt
pip install -e .
