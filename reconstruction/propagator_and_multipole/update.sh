#!/bin/bash
#Usage: ./update.sh [<branch or tag (defaults to 'master')>]
#script credit: https://www.sourcefield.nl/post/git-subtree-survival-tips/
#
REF=${1-master} # branch or tag; defaults to 'master' if parameter 1 not present
REMOTE=remote_repo # just a name to identify the remote
REPO=https://github.com/grantmerz/DESI_Recon.git # replace this with your repository URL
FOLDER=reconstruction/propagator_and_multipole/DESI_Recon # where to mount the subtree

git remote add $REMOTE --no-tags $REPO
if [[ -d $FOLDER ]]; then # update the existing subtree
    git subtree pull $REMOTE $REF --prefix=$FOLDER --squash -m "Merging '$REF' into '$FOLDER'"
else # add the subtree
    git subtree add  $REMOTE $REF --prefix=$FOLDER --squash -m "Merging '$REF' into '$FOLDER'"
fi
git remote remove $REMOTE
