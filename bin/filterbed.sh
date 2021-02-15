#!/bin/bash
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 15/02/2021
# Last Modified Date: 15/02/2021
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>
usage() {
    cat <<HELP_USAGE
    Usage:
    $(basename $0): Filter regions only occurs in n bed files.
    -n [integer] keep regions at least occurs in >= n bed files.
    -o [string] Output file path.
    -b [string] bed file1
    -b [string] bed file2
    ...
    -b [string] bed filen
HELP_USAGE
    exit 0
}
while getopts ":n:o:b:" opt; do
    case $opt in
        n)
            n=${OPTARG}
            ;;
        o)
            output=${OPTARG}
            ;;
        b)
            bedfiles+=("$OPTARG")
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${n}" ] || [ -z "${output}" ] || [ -z "${bedfiles}" ]; then
    usage
    exit
fi

echo "Merging all bed files and keeping regions that occur at least $n files:"
for val in "${bedfiles[@]}"; do
    echo " - $val"
done

cat ${bedfiles[@]} | bedtools sort -i - | bedtools merge -i -  -c 1 -o count | awk -v n=$n '{if($4 >= n ) print}' > $output

