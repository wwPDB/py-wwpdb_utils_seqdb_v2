# OneDep Sequence Search/Retrieval Library

Master: [![Build Status](https://dev.azure.com/wwPDB/wwPDB%20Python%20Projects/_apis/build/status/wwPDB.py-wwpdb_utils_seqdb_v2?branchName=master)](https://dev.azure.com/wwPDB/wwPDB%20Python%20Projects/_build/latest?definitionId=10&branchName=master)

Develop: [![Build Status](https://dev.azure.com/wwPDB/wwPDB%20Python%20Projects/_apis/build/status/wwPDB.py-wwpdb_utils_seqdb_v2?branchName=develop)](https://dev.azure.com/wwPDB/wwPDB%20Python%20Projects/_build/latest?definitionId=10&branchName=develop)

## Introduction

This repository povides the tools and access to remote sequence search and retrieval.

### Installation

Download the library source software from the project repository:

```bash

git clone --recurse-submodules https://github.com/wwpdb/py-wwpdb_utils_seqdb_v2.git

```

Optionally, run test suite using the Tox test runner

```bash
python setup.py test

or simply run

tox

If the environment variable MOCKREQUESTS is set, the tests will not contact any remote servers.


Installation is via the program [pip](https://pypi.python.org/pypi/pip).

```bash
pip install wwpdb.utils.seqdb_v2

or from the local repository:

pip install .
```


