help:
	@echo "Possible targets:"
	@echo "  init: install dependencies using conda or pip"
	@echo "  test: run tests (requires \`nose\`)"
	@echo "  testcov: run tests and determine coverage (requires \`coverage\')"
	@echo "  cov: full coverage report (requires \`coverage\`)"
	@echo "  clean: delete \`.pyc\` files"

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

