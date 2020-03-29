#!/usr/bin/env bash

# 3 opcoes:
#    - $HOME
#    - Caminho completo
#    - apenas com o nome da pas
#    - nomes de ficheiros com "." .
#    - colocar vario fichieros no bash como argumento
#    - (opcao i) interacao


sort_folder() {
    echo "\nEntraste Em: $1 \n"
    cd "$1"
    for FILENAME in * ; do
    
        if [ -f "$FILENAME" ]; then  # if filename is a file, do...
 
            base=$(basename $(pwd))  # current folder name
            name="${FILENAME##*.}"   # file type
            
 
             if [ "$base" = "$name" ]; then
                 continue
             else
                 if [ -d "${FILENAME#*.}" ]; then       # if the file type exists as a folder
                     mv *.${FILENAME#*.} "${FILENAME#*.}"     # Move file to folder with the file type name

                 else
                     # create folder 
                     mkdir "${FILENAME#*.}"      
                     mv *.${FILENAME#*.} "${FILENAME#*.}"

                 fi
             fi
        else
            folder "$FILENAME"
            cd ..
            
        fi
    done
}

sort_folder "$1"



