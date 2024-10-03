#!/bin/bash

# Mostrar todos los archivos .fastq en el directorio actual
echo "Archivos .fastq en el directorio actual:"
ls *.fastq

# Contar el número de líneas en cada archivo .fastq
for file in *.fastq
do
    echo "Contando líneas en $file..."
    wc -l "$file"
done


# Mensaje final
echo "Terminado"
