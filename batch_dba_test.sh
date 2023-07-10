for V in {0..61}
do
  ./dba_test tests/dlbtest_$V.txt tests/dlbtest_$V.log tests/dlbtest_$V.svg
done