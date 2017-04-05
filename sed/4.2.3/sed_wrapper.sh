#!/bin/bash

usage() {
      echo ""
      echo "Usage : sh $0 -e <expression> -i <input file> -o <output file> [-p <posix> -x <extended regex>]"
      echo ""

cat <<'EOF'

  -e <expression>  

  -i </path/to/input file>

  -o </path/to/output file>

  -p posix

  -r extended regular expression

EOF
    exit 0
}

posix=0
ext=0

while getopts ":i:ho:pxe:" opt; do
  case $opt in
    e)
     exp=$OPTARG    
      ;;
    i)
     input=$OPTARG
      ;;
    o)
     output=$OPTARG
      ;;
    p)
     posix=$OPTARG
      ;;
    x)
     ext=$OPTARG
      ;;
    h)
     usage
     exit 1
      ;;    
  esac
done


if [ ! -z "$exp" ] && [ ! -z "$input" ] && [ ! -z "$output" ] && [ "$posix" == 0 ] && [ "$ext" == 0 ]; then
   sed "$exp" "$input" > "$output"

elif [ ! -z "$exp" ] && [ ! -z "$input" ] && [ ! -z "$output" ] && [ "$posix" != 0 ] && [ "$ext" == 0 ]; then
   sed "$exp" "$input" --posix > "$output"

elif [ ! -z "$exp" ] && [ ! -z "$input" ] && [ ! -z "$output" ] && [ "$ext" != 0 ] && [ "$posix" == 0 ]; then
   sed "$exp" "$input" -r > "$output"

elif [ ! -z "$exp" ] && [ ! -z "$input" ] && [ ! -z "$output" ] && [ "$posix" != 0 ] && [ "$ext" != 0 ]; then
   sed "$exp" "$input" --posix -r > "$output"

fi
