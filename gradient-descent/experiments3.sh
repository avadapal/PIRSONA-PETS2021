#./rand_queries 4 &
#./gen_queries &
#./gen_items 4 &
#./gen_users 4 &


#!/bin/bash
for i in `seq 1 10`;
do
        echo $i
        ./p3 localhost > log3;
done    
