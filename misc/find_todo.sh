#!/bin/bash

### ========== DESCRIPTION ========== ###
### Colelcts results, edit path of find and iname 
### ================================= ###

find ./beta*/b*/lDGA_julia -mindepth 1 -maxdepth 3 -iname  'res_NewCode*' > tmp.txt
cat tmp.txt | cut -d / -f 3 > done.txt
find ./beta*/b*/lDGA_julia -mindepth 1 -maxdepth 3 -iname  'config.toml' > tmp.txt
cat tmp.txt | cut -d / -f 3 > all.txt
grep -Fxv -f done.txt all.txt > todo.txt
rm tmp.txt
