#!/bin/bash

# 检查是否提供了文件名参数
if [ -z "$1" ]; then
    echo "Please provide the name of the file containing GEO IDs."
    exit 1
fi

# 读取包含GEO ID的文本文件
while IFS= read -r geo_id
do
    # 检查是否有输入GSE ID号
    if [ -z "$geo_id" ]; then
        echo "Skipping empty line."
        continue
    fi

    # 在当前目录下创建一个文件夹，以GEO ID号命名
    mkdir -p "$geo_id"

    # 递归下载该GSE下的所有元数据
    wget -r -nH --cut-dirs=3 "ftp://ftp.ncbi.nlm.nih.gov/geo/series/${geo_id:0:-3}nnn/$geo_id/"

    # 解压缩下载的压缩文件并删除原始压缩文件
    find . -name "*.zip" -o -name "*.zg" -o -name "*.tgz" -o -name "*.tar" -exec sh -c 'dir="$1"; file="$2"; mkdir -p "${dir%/*}"; tar -xf "$file" -C "${dir%/*}"; rm "$file"' {} "$geo_id" \;

    # 输出下载成功
    echo "Successed downloading and extracting $geo_id"

 # 递归解压缩子文件夹中的 .zip、.zg、.tg、.tar 文件   
function extract_in_subdirectories {
    for folder in "$1"/*; do
        if [ -d "$folder" ]; then
            extract_in_subdirectories "$folder"
        fi
        if [ -f "$folder" ]; then
            if [[ "$folder" == *.zip || "$folder" == *.zg || "$folder" == *.tgz || "$folder" == *.tar ]]; then
                echo "Extracting $folder"
                tar -xf "$folder" -C "$(dirname "$folder")"
                if [ $? -eq 0 ]; then
                    rm "$folder"
                else
                    echo "Failed to extract $folder"
                fi
            fi
        fi
    done
}


    extract_in_subdirectories "$geo_id"

done < "$1"


