#!/bin/sh
git push origin :refs/tags/cpponly
git tag -fa cpponly
git push origin master --tags
