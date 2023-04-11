git add -A
git commit -m "$1"
git tag -d $1
git tag $1
rm -r dist
python -m build
twine upload dist/*
