#!/bin/bash
#Usage: ./update.sh [<branch or tag (defaults to 'master')>]
#script credit: https://www.sourcefield.nl/post/git-subtree-survival-tips/
#
REF=${1-master} # branch or tag; defaults to 'master' if parameter 1 not present
REMOTE=Ginevra_CF # just a name to identify the remote
REPO=https://github.com/gfavole/multi2pcf.git # replace this with your repository URL
FOLDER=clustering_codes/Ginevra_CF/multi2pcf # where to mount the subtree

git remote add $REMOTE --no-tags $REPO
if [[ -d $FOLDER ]]; then # update the existing subtree
    git subtree pull $REMOTE $REF --prefix=$FOLDER --squash -m "Merging '$REF' into '$FOLDER'"
else # add the subtree
    git subtree add  $REMOTE $REF --prefix=$FOLDER --squash -m "Merging '$REF' into '$FOLDER'"
fi
git remote remove $REMOTE
