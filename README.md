# MIDAS (precursor to GAMBIT)

[![Build Status](https://github.com/jlumpe/midas/actions/workflows/ci.yml/badge.svg)](https://github.com/jlumpe/midas/actions/workflows/ci.yml)

MIDAS is an early version of [GAMBIT](https://github.com/jlumpe/gambit), a tool for taxonomic identifiction of bacterial genomes. The final v2.4.0 release was reclassified as GAMBIT v0.1.0 and moved to a new repo. This repository has been kept for reference purposes, the code is way out of date and there's no reason to use it instead of the most recent GAMBIT release.


## Installation

Prerequisites:

    pip install Cython

Clone:

    git clone https://github.com/jlumpe/midas

Build and install:

    cd midas
    python setup.py build_ext --inplace
    python setup.py install
