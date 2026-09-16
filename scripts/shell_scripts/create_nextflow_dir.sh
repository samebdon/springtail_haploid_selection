#!/usr/bin/env bash

NAME="$1"
NF_DIR="nextflow/$NAME"

mkdir -p $NF_DIR/assets
mkdir $NF_DIR/bin
mkdir $NF_DIR/conf
mkdir $NF_DIR/docs
mkdir $NF_DIR/lib
mkdir -p $NF_DIR/modules/local
mkdir $NF_DIR/workflows
touch $NF_DIR/README.md
touch $NF_DIR/main.nf 
touch $NF_DIR/nextflow.config