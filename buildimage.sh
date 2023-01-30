#!/bin/bash

# Build a set of input and output paths for the array_sbatch.sh scripts

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# Copyright 2020-2023, Ben Cardoen


PROJECT="SmlmTools.jl"
set -euo pipefail
NOW=$(date +"%m--%d--%Y ~ %I:%M:%S")
echo "Starting processing at $NOW"

CUR="/home/bcardoen/SFUVault/repositories/$PROJECT"
TMP=/dev/shm
cd $TMP
git clone git@github.com:bencardoen/$PROJECT.git
echo "Creating archive"
zip -rq $PROJECT.zip $PROJECT
rm -rf $TMP/$PROJECT
mv $TMP/$PROJECT.zip $CUR
echo "Done"
cd $CUR

echo "Building Singularity image"
sudo singularity build --force SmlmTools.sif singularity/recipe.def
echo "Done"
