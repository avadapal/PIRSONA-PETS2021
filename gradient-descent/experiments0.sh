for i in `seq 1 10`;
do
        echo $i
        ./p0 localhost localhost localhost $i  > log0;
        sleep 10;
done    

 