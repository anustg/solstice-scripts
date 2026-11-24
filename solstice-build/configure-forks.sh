#!/bin/sh
# Run in the pinned Star-Build checkout. Keep the release dependency graph and
# build flags while selecting public fork URLs and immutable source revisions.
set -eu

use_fork() {
    recipe=$1
    project=$2
    commit=$3
    sed -i -e "s|^tag=.*|tag=\"$commit\"|" \
           -e "s|^  url=.*|  url=\"https://gitlab.com/jdpipe/$project.git\"|" "$recipe"
}

# quadricxy includes the target-aligned parser and quadric geometry support.
use_fork "src/rad-apps/solstice_${SOLSTICE_VERSION}.sh" solstice f5ee7a794421527710775d0c7c16ff1baf3f157e
use_fork src/rad-libs/sanim_0.4.sh solstice-anim 78ac19fc612e7f5828f874423c85e6ba2ca5fcbf
use_fork src/rad-libs/ssolver_0.11.sh solstice-solver 24639a4358e6d770d0a87b9828b0267bb2f47af1
