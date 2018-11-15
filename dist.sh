source /home/u5376227/.virtualenvs/paramagpyenv/bin/activate
python setup.py sdist bdist_wheel
twine upload dist/*