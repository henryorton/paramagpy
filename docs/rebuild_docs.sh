source /home/u5376227/.virtualenvs/paramagpyenv/bin/activate
/home/u5376227/Dropbox/PhD/git/paramagpy/clean.sh
rm -rf build
rm -r ./source/reference/generated/
make html
make latex
