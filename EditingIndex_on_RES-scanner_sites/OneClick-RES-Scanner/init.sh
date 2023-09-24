#!/usr/bin/env bash

# if  argparse not in home directory: - clone it
if [ ! -d ~/argparse-bash ]; then 
 git clone https://github.com/nhoffman/argparse-bash.git $HOME/argparse-bash
fi
git clone https://github.com/ZhangLabSZ/RES-Scanner.git 

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
