#!/bin/bash

usage() {
      echo ""
      echo "Usage : sh $0 -u username -p password -i id -o output"
      echo ""

cat <<'EOF'
  -u username
  -p password
  -i synapse id (one) or file containing multiple id's
  -o output
  -h Show this usage information
EOF
    exit 0
}

while getopts ":u:hp:i:o:" opt; do
  case $opt in
    u)
     username=$OPTARG
      ;;
    p)
     password=$OPTARG
      ;;
    i)
     id=$OPTARG
      ;;
    h)
     usage
     exit 1
      ;;    
    o)
     output=$OPTARG
      ;;  
    \?)
     echo "Invalid option: -$OPTARG" >&2
     exit 1
      ;;
    :)
     echo "Option -$OPTARG requires an argument." >&2
     exit 1
      ;;
  esac
done

if [[ -f $id ]]; then

  while IFS= read -r line

    do 

      synapse -u $username -p $password get $line --downloadLocation $output

    done < $id

else

  synapse -u $username -p $password get $id --downloadLocation $output

fi