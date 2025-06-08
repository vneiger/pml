#!/usr/bin/env bash
#
## 
# this file is adapted from FLINT's make_dist.sh (June 2025)
##
# This script is called by `make dist` 
# in order to create a pml release tarball. It must be called from within
# the root directory of the pml source tree.
#
set -ex

# first argument: the PML version in the form 1.2.3 or 1.2.3-something
# if not given, uses the version in the VERSION file
pml_version=$1

# second, optional argument: the git revision from which to make the
# release; default is to use the HEAD commit.
git_ref=${2:-HEAD}

# prefix used for the content of the tarball, and also the basename of
# the final archive.
archive_prefix="pml-$pml_version"

echo "Exporting from git"
git archive --format tar.gz --prefix "${archive_prefix}/" ${git_ref} > ${archive_prefix}.tar.gz

echo "Extracting"
tar -xf ${archive_prefix}.tar.gz
rm ${archive_prefix}.tar.gz

# update VERSION file
printf $pml_version > VERSION

echo "Bootstrapping"
./bootstrap.sh

echo "Adding / patching / removing files"
# copy some files that should be included in the distribution archive
cp -r config ${archive_prefix}/
cp configure ${archive_prefix}/
cp src/config.h.in ${archive_prefix}/src/
cp VERSION ${archive_prefix}/

# remove some things we don't want to install
pushd ${archive_prefix}
rm -rf .[a-z]*  # no dot files
rm -rf dev

# return to top directory
popd

# create the source archives
echo "Create .tar.gz"
tar -cvzf ${archive_prefix}.tar.gz ${archive_prefix}

echo "Create .tar.xz"
tar -cJf ${archive_prefix}.tar.xz ${archive_prefix}

echo "Create .zip"
zip -9 -r ${archive_prefix}.zip ${archive_prefix}
