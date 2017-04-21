# GeCKOReader.jl

![Project Status](https://img.shields.io/badge/status-wip-red.svg)[![Build
Status](https://travis-ci.org/tlnagy/GeCKOReader.jl.svg?branch=master)](https://travis-ci.org/tlnagy/GeCKOReader.jl)

Some analysis code for interpreting GeCKOv2 CRISPR knockout screens

## Loading

These functions are useful for loading `fastq.gz` files

```@autodocs
Modules=[GeCKOReader]
Private=false
Pages=["io.jl", "utils.jl"]
```

## Processing

```@autodocs
Modules=[GeCKOReader]
Private=false
Pages=["processing.jl"]
```

## Plotting

```@autodocs
Modules=[GeCKOReader.Plot]
Pages=["plotting.jl"]
```
