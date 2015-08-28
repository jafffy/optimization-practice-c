for ((i = 0; i < 10; ++i)) do
START=$(date +%s.%N)
./serial
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo $DIFF >> benchmark-opt.txt
done
