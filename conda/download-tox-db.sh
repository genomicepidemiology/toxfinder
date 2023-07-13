#!/usr/bin/env bash

echo "Downloading lastest version of the toxfinder database to current directory..."

mkdir toxfinder_db
cd toxfinder_db

wget https://bitbucket.org/genomicepidemiology/toxfinder_db/get/master.tar.gz
tar -xvf master.tar.gz --strip-components 1

echo "Installing the ToxFinder database with KMA"
python INSTALL.py

echo "The ToxFinder database has been downloaded and installed."

exit 0
