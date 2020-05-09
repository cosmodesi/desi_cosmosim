#!/bin/bash
#Usage: ./update.sh [<branch or tag (defaults to 'master')>]
#script credit: https://www.sourcefield.nl/post/git-subtree-survival-tips/
#
REF=${1-master} # branch or tag; defaults to 'master' if parameter 1 not present
REMOTE=ChengZhao # just a name to identify the remote
REPO=https://github.com/cheng-zhao/jobfork.git # replace this with your repository URL
FOLDER=massive_production/MPI_wrapper/jobfork # where to mount the subtree

git remote add $REMOTE --no-tags $REPO
if [[ -d $FOLDER ]]; then # update the existing subtree
    git subtree pull $REMOTE $REF --prefix=$FOLDER --squash -m "Merging '$REF' into '$FOLDER'"
else # add the subtree
    git subtree add  $REMOTE $REF --prefix=$FOLDER --squash -m "Merging '$REF' into '$FOLDER'"
fi
git remote remove $REMOTE
