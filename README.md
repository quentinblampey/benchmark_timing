# Timing three public algorithms

This repository contains code from ACDC, MP, and Phenograph, three public algorithms. The purpose is to time these algorithms on different dataset sizes.

## Getting started

Create a fresh Python 3.7 environment and install the dependencies:

```
pip install -r requirements.txt
```

## Run the scripts (on CPUs)

Run Phenograph (expected: about one day)

```
python timing_pheno.py
```

Run ACDC (expected: about one day)

```
python timing_acdc.py
```

Run ACDC2 (expected: about one day)

```
python timing_acdc2.py
```

Run MP (expected: about one day)

```
python timing_mp.py
```
