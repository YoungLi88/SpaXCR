.PHONY: docs docs-clean

docs:
	sphinx-build -b html docs docs/_build/html

docs-clean:
	rm -rf docs/_build
