build:
	make test && make clean && tar -cvzf pymolprobity.tar.gz --exclude=.DS_Store pymolprobity

init:
	conda env create -f environment.yml || pip install -r requirements.txt

test:
	export PYTHONPATH=${PYMOL_HOME}; nosetests tests

testcov:
	export PYTHONPATH=${PYMOL_HOME}; nosetests --with-coverage --cover-erase tests

cov:
	make testcov
	coverage report -m pymolprobity/*.py

clean:
	rm pymolprobity/*.pyc
	rm tests/*.pyc

