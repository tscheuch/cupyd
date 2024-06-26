POETRY_VERSION = 1.2.2
PYTHON_VERSION = 3.9

# Managing
# ========
.PHONY: poetry
poetry:
	# Please make sure the <python> interpreter is pointing to Python 3.7 or higher.
	# If that's not the case, run this command manually by setting a proper version.
	curl -sSL https://install.python-poetry.org | python - --version $(POETRY_VERSION)

.PHONY: venv-with-dependencies
venv-with-dependencies:
	python$(PYTHON_VERSION) -m venv .venv
	poetry run pip install --upgrade pip
	poetry install

# Linting
# =======
.PHONY: black
black:
	poetry run black --check .

.PHONY: black!
black!:
	poetry run black .

.PHONY: isort
isort:
	poetry run isort . --check --diff

.PHONY: isort!
isort!:
	poetry run isort .

.PHONY: mypy
mypy:
	poetry run mypy cupyd/georef.py

# Testing
# =======
.PHONY: test-results
test-results:
	python test_llanquihue_results.py
