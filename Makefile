POETRY_VERSION = 1.2.0b2
PYTHON_VERSION = 3.8

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

# GHA Workflows
.PHONY: build
build:
	poetry install --no-dev

.PHONY: build-dev
build-dev:
	poetry install


# LINTERS
black:
	poetry run black swmm_modflow_oss --check

black!:
	poetry run black swmm_modflow_oss

isort:
	poetry run isort swmm_modflow_oss/* --check --settings-path ./pyproject.toml --diff

isort!:
	poetry run isort swmm_modflow_oss/* --settings-path ./pyproject.toml
