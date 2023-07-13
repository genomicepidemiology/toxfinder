#!/bin/bash

mkdir -p ${PREFIX}/bin

chmod +x toxfinder-1.0.py
cp toxfinder-1.0.py ${PREFIX}/bin/toxfinder-1.0.py

# copy script to download database
chmod +x ${RECIPE_DIR}/download-tox-db.sh
cp ${RECIPE_DIR}/download-tox-db.sh ${PREFIX}/bin/download-tox-db.sh
