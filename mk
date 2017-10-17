ofiles=`echo *.o`

rm *.o

make -f Makefile > make.out 2> make.err

rm *.o

cat make.err
ls -l muscle
