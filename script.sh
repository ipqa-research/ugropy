#!/bin/bash

# Define la ruta del directorio que contiene los archivos .py
dir="./tests/dortmund"

# Itera sobre todos los archivos .py en el directorio
for file in "$dir"/*.py; do
    # Realiza el reemplazo en cada archivo utilizando sed
    sed -i 's/assert ug.get_dortmund_groups(identifier, identifier_type) == result/assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result/g' "$file"
    # sed -i 's/assert ug.get_psrk_groups(identifier, identifier_type) == result/assert ug.get_groups(ug.psrk, identifier, identifier_type) == result/g' "$file"
done