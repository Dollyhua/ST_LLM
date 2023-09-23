#!/bin/bash

# 判断是否有输入GSE ID号
if [ -z "$1" ]; then
    echo "Please provide a GSE number."
    exit 1
fi

# 在当前目录下创建一个文件夹，以GSE ID号命名
mkdir -p "$1"

# 递归下载该GSE下的所有元数据
wget -r -nH --cut-dirs=3 "ftp://ftp.ncbi.nlm.nih.gov/geo/series/${1:0:-3}nnn/$1/"

# 解压缩下载的压缩文件并删除原始压缩文件
find . -name "*.zip" -o -name "*.zg" -o -name "*.tgz" -o -name "*.tar" -exec sh -c 'dir="$1"; file="$2"; mkdir -p "${dir%/*}"; tar -xf "$file" -C "${dir%/*}"; rm "$file"' {} "$1" \;

# 输出下载成功
echo "Successed downloading and extracting $1"

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
                rm "$folder"
            fi
        fi
    done
}

extract_in_subdirectories "$1"