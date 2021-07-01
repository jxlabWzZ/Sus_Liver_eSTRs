for i in `cat chip7.id`
do
g=`echo $i | sed 's/_/\t/g' | cut -f 1`
c=`echo $i | sed 's/_/\t/g' | cut -f 2`
p=`echo $i | sed 's/_/\t/g' | cut -f 3`

sed "s/ENSSSCG00000003653/$g/g;s/chr6/$c/g;s/94905102/$p/g;s/plot\/f7/plot\/chipseq77\/f7/g" dd.py > chip77.py
chmod 755 chip77.py

nohup python chip77.py &

sleep 15
done
