#!/bin/bash

cd singularity_img

FILES=$(ls *.def)

for i in $FILES; do
	IMG_FN="${i%.def}.simg"
	if [ -f $IMG_FN ]; then
		echo "$i already present, please, delete first if you want to re-build"
	else
		echo "building $i"
		sudo singularity build $IMG_FN $i
	fi
done

cd ..
