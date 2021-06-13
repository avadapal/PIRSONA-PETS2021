#./rand_queries 4 &
#./gen_queries &
#./gen_items 4 &
#./gen_users 4 &
 
for i in `seq 1 10`;
do
        echo $i
        ./p1 localhost localhost localhost $i > log1;
        sleep 4;
done   
