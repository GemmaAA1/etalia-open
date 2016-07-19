### Install/Update
Check out `./scripts/routines/update.py -h` for help. Examples use cases:

* `./scripts/routines/update.py -i`: (init) Populate database and initialize NLP Engines
* `./scripts/routines/update.py -l [file_name (optional)]`: (load) Load data from fixture

### After rebase/merge :

* Update your common.py file based on diff from common.py.dist
* Run `./scripts/routines/update.py`
