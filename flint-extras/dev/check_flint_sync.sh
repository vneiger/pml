#!/bin/sh
#
# Copyright (C) 2026 Vincent Neiger
#
# This file is part of PML.
#
# PML is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
# as published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version. See
# <https://www.gnu.org/licenses/>.
#
# Compares the parts of PML's build system that are copies of FLINT's with a
# FLINT checkout, see dev/SYNC_WITH_FLINT.md.
#
# Usage, from flint-extras/:
#     dev/check_flint_sync.sh [-v] FLINT_DIR
#
# Prints one line per part, "ok" or "DRIFT"; with -v, also the diff of each
# part that drifted (PML's version first). The exit status is 1 if some part
# drifted, 2 on a usage error.

verbose=0
if test "$1" = "-v"; then verbose=1; shift; fi

if test $# -ne 1 || ! test -f "$1/configure.ac" || ! test -f configure.ac;
then
    echo "usage: dev/check_flint_sync.sh [-v] FLINT_DIR  (run from flint-extras/)" >&2
    exit 2
fi

flint=$1
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
status=0

# report NAME PML_FILE FLINT_FILE
report ()
{
    if cmp -s "$2" "$3";
    then
        printf '  ok     %s\n' "$1"
    else
        printf '  DRIFT  %s\n' "$1"
        status=1
        if test $verbose -eq 1; then diff -u "$2" "$3" | sed 's/^/         /'; fi
    fi
}

# section SECTION FILE: the body of the section of the configure.ac FILE
# whose banner title starts with SECTION, up to the banner of the next
# section; the issue URL is normalised so that PML's and FLINT's versions
# compare equal
section ()
{
    awk -v title="# $1" 'index($0, title) == 1 { s = 1; next }
         s == 1 && /^####/ { s = 2; next }
         s == 2 && /^####/ { exit }
         s == 2 { print }' "$2" |
    sed 's|https://github.com/vneiger/pml/issues/|ISSUES_URL|;
         s|https://github.com/flintlib/flint/issues/|ISSUES_URL|'
}

# the definitions of the host patterns (X86_64_PATTERN, ...) of acinclude.m4
patterns ()
{
    grep -A1 '^define([A-Z0-9_]*_PATTERN,' "$1" | grep -v '^--$'
}

echo "Comparing with FLINT at $flint ($(git -C "$flint" log -1 --format='%h %cs' 2>/dev/null || echo 'not a git checkout'))"

for f in config.guess config.sub configfsf.guess configfsf.sub;
do
    report "config/$f" "config/$f" "$flint/config/$f"
done

patterns acinclude.m4 > "$tmp/pml_patterns"
patterns "$flint/acinclude.m4" > "$tmp/flint_patterns"
report "acinclude.m4: host patterns" "$tmp/pml_patterns" "$tmp/flint_patterns"

for sec in "architecture specifics" "CFLAGS";
do
    section "$sec" configure.ac > "$tmp/pml_section"
    section "$sec" "$flint/configure.ac" > "$tmp/flint_section"
    report "configure.ac: section $sec" "$tmp/pml_section" "$tmp/flint_section"
done

exit $status
