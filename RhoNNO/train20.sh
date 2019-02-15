for i in 20 21 22 23 24 25 26 27; do
cat <<EOF > trackml$i.nno
xmlp 9 15 5 1
transfer TR_FERMI
momentum 0.2
#balance true
plots false
test 10000
start 1
stop 350
target 5.0
tree tracks$i
#input 0.001*rz1:phi1:0.001*z1:0.001*rz2:phi2:0.001*z2:f0:f1:0.001*f2:0.001*f3:0.001*f4:0.001*f5
input 0.001*rz1:phi1:0.001*z1:0.001*rz2:phi2:0.001*z2:f0:f1:0.001*score
#input f0:f1:0.001*f2:0.001*f3:0.001*f4:0.001*f5
output truth==2

#define the data source
datapath /Users/marcel/workspace/training_xmlp
networkpath ../Networks$i

#file event000021100.root
file event000021101.root
file event000021102.root
file event000021103.root
file event000021104.root
file event000021105.root
file event000021106.root
file event000021107.root
file event000021108.root
file event000021109.root
file event000021110.root
file event000021111.root
file event000021112.root
file event000021113.root
file event000021114.root
file event000021115.root
file event000021116.root
file event000021117.root
file event000021118.root
file event000021119.root
file event000021120.root
file event000021121.root
file event000021122.root
file event000021123.root
file event000021124.root
file event000021125.root
file event000021126.root
file event000021127.root
file event000021128.root
file event000021129.root
file event000021130.root
file event000021131.root
file event000021132.root
file event000021133.root
file event000021134.root
file event000021135.root
file event000021136.root
file event000021137.root
file event000021138.root
file event000021139.root
file event000021140.root
file event000021141.root
file event000021142.root
file event000021143.root
file event000021144.root
file event000021145.root
file event000021146.root
file event000021147.root
file event000021148.root
file event000021149.root
file event000021150.root
file event000021151.root
file event000021152.root
file event000021153.root
file event000021154.root
file event000021155.root
file event000021156.root
file event000021157.root
file event000021158.root
file event000021159.root
file event000021160.root
file event000021161.root
file event000021162.root
file event000021163.root
file event000021164.root
file event000021165.root
file event000021166.root
file event000021167.root
file event000021168.root
file event000021169.root
EOF

../bin/NetworkTrainer trackml$i.nno >trackml$i.out
cp ../Networks$i/NNO0000.TXMLP XMLP$i.net

tail trackml$i.out | mail -s trackml$i.out marcel@cloudkitchen.info

done


