# Develop Guidelines
## environment prepare
```
conda create -n trackc python=3.9
conda activate trackc
pip install poetry
poetry install
```
## pre commit
```
# Format python code
poe format
# Run test code
poe test
# Local running document
poe docs-serve
```
## release and pubish
```
git tag vX.X.X
git push origin vX.X.X
```