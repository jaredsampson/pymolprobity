build:
	make test && make clean && tar -cvzf pymolprobity.tar.gz --exclude=.DS_Store pymolprobity

init:
	conda env create -f environment.yml || pip install -r requirements.txt

updateenv:
	conda env export -f environment.yml
	which -s pip && pip freeze > requirements.txt

test:
	$(eval export PYTHONPATH=${PYMOL_HOME})
	nosetests tests

testcov:
	$(eval export PYTHONPATH=${PYMOL_HOME})
	nosetests --with-coverage --cover-erase

cov:
	make testcov
	coverage report -m pymolprobity/*.py

clean:
	rm pymolprobity/*.pyc
	rm tests/*.pyc

run:
ifndef pml
	echo "No test indicated, just starting PyMOL and importing pymolprobity."
	pymol pml/test0.pml
else
	pymol pml/test${pml}.pml
endif
