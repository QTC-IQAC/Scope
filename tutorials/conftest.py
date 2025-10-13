# conftest.py
import os
import nbformat

def clear_notebook_outputs(notebook_path):
    nb = nbformat.read(notebook_path, as_version=nbformat.NO_CONVERT)
    for cell in nb.cells:
        if "outputs" in cell:
            cell["outputs"] = []
        if "execution_count" in cell:
            cell["execution_count"] = None
    nbformat.write(nb, notebook_path)

def pytest_sessionstart(session):
    """Hook executed before any tests run"""
    tutorials_dir = os.path.abspath("./tutorials/")  # adjust if needed
    print(tutorials_dir)
    for nb_file in [file for file in os.listdir(tutorials_dir) if file.endswith(".ipynb")]: #and 'tests' not in file]:
        path = f"{tutorials_dir}/{nb_file}"
        clear_notebook_outputs(path)
