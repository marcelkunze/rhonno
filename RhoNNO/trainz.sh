for i in 20 21 22 23 24 25 26 27; do
cat <<EOF > trainz$i.nno
xmlp 9 25 15 1
transfer TR_LINEAR
momentum 0.2
#balance true
plots false
test 10000
start 1
stop 200
target 0.2
tree tracks$i
input 0.001*rz1:phi1:0.001*z1:0.001*rz2:phi2:0.001*z2:f0:f1:0.001*f4:0.001*f3
output 0.001*abs(vz)
cut l1!=l2&&truth!=0

#define the data source
datapath /Users/marcel/workspace/training_octant
networkpath ../Networks$i

#file event000021100.root
file event000021101.root
file event000021102.root
#file event000021103.root
#file event000021104.root
file event000021105.root
file event000021106.root
#file event000021107.root
#file event000021108.root
#file event000021109.root
EOF

../bin/NetworkTrainer trainz$i.nno >trainz$i.out
cp ../Networks$i/NNO0000.TXMLP XMLPV$i.net

tail trainz$i.out | mail -s trainz$i.out marcel@cloudkitchen.info

done


