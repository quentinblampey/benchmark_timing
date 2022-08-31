# Timing three public algorithms

This repository contains code from ACDC, MP, and Phenograph, three public algorithms. The purpose is to time these algorithms on different dataset sizes.

## Getting started

Create two fresh Python 3.7 environments and install dependencies:

- On the first env, run `pip install -r requirements.txt`
- On the second env, run `pip install -r requirements_mp.txt`

## Run the jobs (on CPUs)

Three jobs have to be run: `acdc.sh`, `mp.sh`, and `pheno.sh`.

> Note that, you have to activate the first env, except for `mp.sh` (second env).

All we need are the logs, which tells how long the algorithms run on the different dataset sizes.
