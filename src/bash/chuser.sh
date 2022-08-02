#!/bin/bash

TMP="some_temporary_filename_this_is_dumb"

RECURSIVE=""
POSITIONAL=()
while [[ $# -gt 0 ]]
do
	key="$1"

	case $key in
		-r|--recursive)
			RECURSIVE="-r"
			shift # past argument
			;;
		*)    # unknown option
			POSITIONAL+=("$1") # save it in an array for later
			shift # past argument
			;;
	esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

cp -d --preserve=all $RECURSIVE $1 $TMP || exit 1
rm $RECURSIVE $1
mv $TMP $1
