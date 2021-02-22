#!/bin/bash

find $1 -size -99M -type f -print0 | xargs -0 git add
git commit -m $2
git push origin master
