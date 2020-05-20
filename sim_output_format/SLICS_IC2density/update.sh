#!/bin/bash
#Usage: ./update.sh [<branch or tag (defaults to 'master')>]
#script credit: https://www.sourcefield.nl/post/git-subtree-survival-tips/
#
REF=${1-master} # branch or tag; defaults to 'master' if parameter 1 not present
REMOTE=slics2den # just a name to identify the remote
REPO=https://github.com/balaguera/DESI_balaguera.git # replace this with your repository URL
FOLDER=sim_output_format/SLICS_IC2density/DESI_balaguera # where to mount the subtree

git remote add $REMOTE --no-tags $REPO
if [[ -d $FOLDER ]]; then # update the existing subtree
    git subtree pull $REMOTE $REF --prefix=$FOLDER --squash -m "Merging '$REF' into '$FOLDER'"
else # add the subtree
    git subtree add  $REMOTE $REF --prefix=$FOLDER --squash -m "Merging '$REF' into '$FOLDER'"
fi
git remote remove $REMOTE
