# Develop Guidelines
## environment prepare
```
conda create -n trackc python=3.8
conda activate trackc
pip install poetry
poetry install
```
## pre commit
```
poe format
poe test
```
## release and pubish
```
git tag vX.X.X
git push origin vX.X.X
```