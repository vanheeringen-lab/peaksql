from setuptools import setup
import toml

project = toml.load("pyproject.toml")["project"]

# read the readme as long description
with open("README.md") as f:
    project["long_description"] = f.read()

project["long_description_content_type"] = "text/markdown"

setup(**project)
