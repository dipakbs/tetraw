for i in {1..400}
do
./a.out input_${i}.dat out_{i}.dat
grep "Tetrhedral " out_{i}.dat 
done > order.dat
