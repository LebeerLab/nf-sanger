fin=$1
fout=$2
  
while read line ; do 
echo $line
read line ; echo $line | rev
read line ; echo $line
read line ; echo $line | rev 
done < $fin >$fout
