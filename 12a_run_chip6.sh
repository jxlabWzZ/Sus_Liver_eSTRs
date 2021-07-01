for i in `cat chip6.id`
do
g=`echo $i | sed 's/_/\t/g' | cut -f 1`
c=`echo $i | sed 's/_/\t/g' | cut -f 2`
p=`echo $i | sed 's/_/\t/g' | cut -f 3`

sed "s/ENSSSCG00000003653/$g/g;s/chr6/$c/g;s/94905102/$p/g;s/plot\/f6/plot\/chipseq66\/f6/g" bb.py > chip66.py
chmod 755 chip66.py

nohup python chip66.py &

sleep 15
done
