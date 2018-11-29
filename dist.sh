source /home/u5376227/.virtualenvs/paramagpyenv/bin/activate
rm -r dist
rm -r build
python setup.py sdist bdist_wheel
twine upload dist/*
